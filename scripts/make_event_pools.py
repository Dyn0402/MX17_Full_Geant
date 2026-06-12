#!/usr/bin/env python3
"""
make_event_pools.py — PLAN Stage 3: slim Geant4 output into event pools
========================================================================
Reads HitTree + EventTree from mx17_full_sim output and writes compact npz
pools for the Python pile-up sampler (MX17_Simulation, PLAN Stage 4).

Per-step HitTree rows are aggregated into one DIGEST row per
(event, track-class, arm, detector):

    edep_sum, n_steps, t_first/t_last [ns, Geant4 global], ke_first,
    entry (u,v,w) and exit (u,v,w) positions [mm]

which preserves everything the fast-MC needs (trigger logic, MM entry point +
direction, calorimeter sums, timing) at ~50 B/row instead of ~100 raw steps.

Modes
-----
--mode pairs    : input = X17+IPC pair runs (event_type 0/1).
                  Output: pool_x17.npz, pool_ipc.npz.  Events with >=1 hit.
                  Truth: vertex, KEs, direction unit vectors, opening angle,
                  inv_mass.  Pairs are generated at t=0 (no E_n): the sampler
                  assigns pulse times from the measured capture-rate curve
                  (analysis/mev/mev_rates.json).
--mode neutrons : input = neutron-beam runs (event_type 2).
                  Output: pool_neutron_bg.npz.  Events with >=1 hit.
                  Truth: neutron_E_eV, capture volume/process (coded), cap pos.
                  Pulse time of every hit = TOF(E_n over L_eff) + Geant4 time;
                  TOF is computed here and stored per event (t0_pulse_ns).
                  Optional --t-max-us keeps only events whose FIRST hit falls
                  inside the pulse window of interest (late thermal captures
                  are better modelled as a steady rate than event-by-event).

Normalisation (stored in meta): each pool event represents
    n_per_pulse / n_events_simulated
beam neutrons (neutron mode), so   rate/pulse = n_pool * n_per_pulse / n_sim.
Pairs pools are normalised downstream by the capture rate instead.

Usage (lxplus):
    python3 scripts/make_event_pools.py --mode pairs \
        /eos/experiment/ntof/data/x17/full_sim/pairs_v2_step_target/ \
        -o analysis/pools/ --workers 6
    python3 scripts/make_event_pools.py --mode neutrons \
        /eos/experiment/ntof/data/x17/full_sim/neutrons_fullrange/ \
        -o analysis/pools/ --workers 6 --t-max-us 50
"""

import argparse
import glob
import json
import sys
from datetime import date
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

# ── Constants ─────────────────────────────────────────────────────────────────
L_FLIGHT_M   = 19.5     # EAR2 nominal flight path to target [m]
GUN_Y_M      = 0.20     # neutron gun launches at y = -20 cm (t=0 there)
M_N_EV       = 939.565e6
C_M_S        = 2.998e8
N_PER_PULSE  = 2.2628e7  # flux-file integral 1 meV-100 MeV (fullrange runs)

DET_CODE = {"DriftGas": 0, "AmpGas": 1, "PlasticScint": 2, "LiqScint_1": 3,
            "LiqScint_2": 4, "BackScintL": 5, "BackScintR": 6}
CAP_VOL_CODE = {"": 0, "He3Gas": 1, "He3Cap_Al": 2, "He3Cap_CFRP": 3,
                "PlasticScint": 4, "LiqScint_1": 5, "LiqScint_2": 6,
                "BackScintL": 7, "BackScintR": 8, "World": 9}   # 10 = other
CAP_PROC_CODE = {"": 0, "neutronInelastic": 1, "nCapture": 2}    # 3 = other

DIGEST_COLS = ["evt_idx", "trackClass", "armID", "det", "n_steps",
               "edep_MeV", "t_first_ns", "t_last_ns", "ke_first_MeV",
               "u_first", "v_first", "w_first", "u_last", "v_last", "w_last"]

PAIR_TRUTH = ["event_type", "vtx_x", "vtx_y", "vtx_z", "inv_mass",
              "em_ke", "ep_ke", "em_px", "em_py", "em_pz",
              "ep_px", "ep_py", "ep_pz", "openingAngle"]


def tof_ns(E_eV, L_m):
    """Relativistic neutron time of flight [ns] over L_m metres."""
    gamma = 1.0 + np.asarray(E_eV, dtype=float) / M_N_EV
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    return L_m / (beta * C_M_S) * 1e9


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+",
                   help="ROOT files, globs, or directories")
    p.add_argument("--mode", choices=["pairs", "neutrons"], required=True)
    p.add_argument("-o", "--outdir", default="analysis/pools",
                   help="Output directory for pool npz files")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--max-files", type=int, default=None)
    p.add_argument("--chunk-size", type=int, default=500_000,
                   help="HitTree rows per uproot chunk")
    p.add_argument("--t-max-us", type=float, default=None,
                   help="neutrons mode: keep only events whose first hit has "
                        "pulse time (TOF + Geant4) below this [us]")
    p.add_argument("--max-events-per-class", type=int, default=None,
                   help="pairs mode: random-subsample each class to this size")
    p.add_argument("--flight-path-m", type=float, default=L_FLIGHT_M)
    return p.parse_args()


def resolve_files(inputs, max_files):
    files = []
    for inp in inputs:
        path = Path(inp)
        if path.is_dir():
            files += sorted(glob.glob(str(path / "*.root")))
        else:
            files += sorted(glob.glob(inp))
    files = sorted(set(files))
    return files[:max_files] if max_files else files


# ── Per-file digest (worker) ──────────────────────────────────────────────────
def digest_file(args):
    """Stream one file's HitTree -> per-(event,trackClass,arm,det) digests.

    Returns dict with 'digests' (DataFrame), 'events' (EventTree DataFrame for
    hit-bearing events only), or {'err': ...}.  eventIDs are file-local; the
    merger re-indexes them globally.
    """
    path, mode, chunk_size = args
    try:
        with uproot.open(path) as f:
            evt_branches = (PAIR_TRUTH + ["eventID"] if mode == "pairs" else
                            ["eventID", "neutron_E_eV", "capture_vol",
                             "capture_proc", "cap_y"])
            evts = f["EventTree"].arrays(evt_branches, library="pd")

            hit_branches = ["eventID", "trackID", "armID", "detType",
                            "u", "v", "w", "edep", "ke", "time"]
            parts, carry = [], None
            for chunk in f["HitTree"].iterate(hit_branches,
                                              step_size=chunk_size,
                                              library="pd"):
                if carry is not None:
                    chunk = pd.concat([carry, chunk], ignore_index=True)
                # hold back the (possibly incomplete) last event of the chunk
                last_eid = chunk["eventID"].iloc[-1]
                carry = chunk[chunk["eventID"] == last_eid]
                cur = chunk[chunk["eventID"] != last_eid]
                if len(cur):
                    parts.append(_digest_chunk(cur))
            if carry is not None and len(carry):
                parts.append(_digest_chunk(carry))

        if not parts:
            return {"path": path, "digests": None, "events": None}
        dig = pd.concat(parts, ignore_index=True)
        keep = evts[evts["eventID"].isin(dig["eventID"].unique())].copy()
        print(f"[done] {Path(path).name}: {len(keep):,} hit-bearing events, "
              f"{len(dig):,} digest rows", flush=True)
        return {"path": path, "digests": dig, "events": keep,
                "n_events_total": len(evts)}
    except Exception as e:
        print(f"[SKIP] {Path(path).name}: {e.__class__.__name__}: "
              f"{str(e)[:90]}", flush=True)
        return {"path": path, "err": str(e)[:200]}


def _digest_chunk(df):
    """Aggregate per-step rows -> digest rows (complete events only)."""
    det = df["detType"].map(lambda s: DET_CODE.get(str(s), -1)).astype(np.int8)
    tid = df["trackID"].to_numpy()
    track_class = np.where(tid == 1, 0, np.where(tid == 2, 1, 2)).astype(np.int8)
    g = pd.DataFrame({
        "eventID": df["eventID"].to_numpy(),
        "trackClass": track_class,
        "armID": df["armID"].astype(np.int8),
        "det": det,
        "edep": df["edep"].to_numpy() * 1e-6,        # eV -> MeV
        "ke": df["ke"].to_numpy(),
        "t": df["time"].to_numpy(),
        "u": df["u"].to_numpy(), "v": df["v"].to_numpy(),
        "w": df["w"].to_numpy(),
    })
    g = g[g["det"] >= 0]
    if not len(g):
        return pd.DataFrame()
    # rows are time-ordered within a track's traversal as written; use
    # first/last by original order for entry/exit, min/max for times
    gb = g.groupby(["eventID", "trackClass", "armID", "det"], sort=False)
    out = gb.agg(
        n_steps=("edep", "size"), edep_MeV=("edep", "sum"),
        t_first_ns=("t", "min"), t_last_ns=("t", "max"),
        ke_first_MeV=("ke", "first"),
        u_first=("u", "first"), v_first=("v", "first"), w_first=("w", "first"),
        u_last=("u", "last"), v_last=("v", "last"), w_last=("w", "last"),
    ).reset_index()
    return out


# ── Pool assembly ─────────────────────────────────────────────────────────────
def assemble_pool(results, mode, args, class_filter=None):
    """Merge per-file digests into flat CSR arrays for one pool."""
    ev_frames, dig_frames, n_sim = [], [], 0
    next_idx = 0
    for r in results:
        if r.get("digests") is None:
            continue
        evts, dig = r["events"], r["digests"]
        n_sim += r["n_events_total"]
        if class_filter is not None:
            evts = evts[evts["event_type"] == class_filter]
            if not len(evts):
                continue
            dig = dig[dig["eventID"].isin(evts["eventID"])]
        # file-local eventID -> global pool index
        idx_map = {eid: next_idx + i for i, eid in
                   enumerate(evts["eventID"].to_numpy())}
        next_idx += len(evts)
        dig = dig.copy()
        dig["evt_idx"] = dig["eventID"].map(idx_map)
        ev_frames.append(evts)
        dig_frames.append(dig.drop(columns=["eventID"]))

    if not ev_frames:
        return None
    events = pd.concat(ev_frames, ignore_index=True)
    digests = pd.concat(dig_frames, ignore_index=True)

    # neutrons: pulse-time anchor per event, optional window filter
    if mode == "neutrons":
        L_eff = args.flight_path_m - GUN_Y_M
        events["t0_pulse_ns"] = tof_ns(events["neutron_E_eV"].to_numpy(), L_eff)
        if args.t_max_us is not None:
            tfirst = digests.groupby("evt_idx")["t_first_ns"].min()
            t_pulse = events["t0_pulse_ns"].to_numpy()
            t_pulse = t_pulse + tfirst.reindex(range(len(events))).to_numpy()
            keep = np.nonzero(t_pulse <= args.t_max_us * 1e3)[0]
            events, digests = _subselect(events, digests, keep)

    # pairs: optional subsample
    if mode == "pairs" and args.max_events_per_class and \
            len(events) > args.max_events_per_class:
        rng = np.random.default_rng(42)
        keep = np.sort(rng.choice(len(events), args.max_events_per_class,
                                  replace=False))
        events, digests = _subselect(events, digests, keep)

    # CSR offsets: sort digests by evt_idx, record starts
    digests = digests.sort_values("evt_idx", kind="stable").reset_index(drop=True)
    counts = np.bincount(digests["evt_idx"].to_numpy(), minlength=len(events))
    dig_start = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)

    return events, digests, dig_start, n_sim


def _subselect(events, digests, keep_idx):
    """Keep events at positions keep_idx; remap digest evt_idx."""
    remap = -np.ones(len(events), dtype=np.int64)
    remap[keep_idx] = np.arange(len(keep_idx))
    events = events.iloc[keep_idx].reset_index(drop=True)
    digests = digests.copy()
    digests["evt_idx"] = remap[digests["evt_idx"].to_numpy()]
    digests = digests[digests["evt_idx"] >= 0]
    return events, digests


def save_pool(path, events, digests, dig_start, n_sim, mode, args, files):
    meta = {
        "created": str(date.today()),
        "mode": mode,
        "source": f"{len(files)} files: {files[0]} ...",
        "geometry": ("STEP capsule: He-3 500 bar (r=10 mm bore, 60 mm on-axis), "
                     "Al 0.6 mm barrel / 5 mm dome + CFRP 0.9 mm, "
                     "MM front face 250 mm, ArIso drift gas"),
        "n_events_simulated": int(n_sim),
        "n_events_pool": int(len(events)),
        "det_code": DET_CODE,
        "cap_vol_code": CAP_VOL_CODE if mode == "neutrons" else None,
        "cap_proc_code": CAP_PROC_CODE if mode == "neutrons" else None,
        "track_class": {"0": "primary e- (trackID 1)",
                        "1": "primary e+ (trackID 2)", "2": "secondary"},
        "units": {"edep": "MeV", "ke": "MeV", "t": "ns", "uvw": "mm"},
        "digest_note": ("one row per (event, trackClass, arm, det); "
                        "t is Geant4 global time (t=0 at primary launch)"),
    }
    arrays = {
        "dig_start": dig_start,
        "dig_trackClass": digests["trackClass"].to_numpy(np.int8),
        "dig_armID": digests["armID"].to_numpy(np.int8),
        "dig_det": digests["det"].to_numpy(np.int8),
        "dig_n_steps": digests["n_steps"].to_numpy(np.uint16),
        "dig_edep_MeV": digests["edep_MeV"].to_numpy(np.float32),
        "dig_t_first_ns": digests["t_first_ns"].to_numpy(np.float32),
        "dig_t_last_ns": digests["t_last_ns"].to_numpy(np.float32),
        "dig_ke_first_MeV": digests["ke_first_MeV"].to_numpy(np.float32),
        "dig_u_first": digests["u_first"].to_numpy(np.float32),
        "dig_v_first": digests["v_first"].to_numpy(np.float32),
        "dig_w_first": digests["w_first"].to_numpy(np.float32),
        "dig_u_last": digests["u_last"].to_numpy(np.float32),
        "dig_v_last": digests["v_last"].to_numpy(np.float32),
        "dig_w_last": digests["w_last"].to_numpy(np.float32),
    }
    if mode == "pairs":
        for col in PAIR_TRUTH:
            arrays[f"evt_{col}"] = events[col].to_numpy(np.float32)
    else:
        meta["n_per_pulse"] = N_PER_PULSE
        meta["rate_per_pulse"] = float(len(events) / n_sim * N_PER_PULSE)
        meta["flight_path_m"] = args.flight_path_m
        meta["gun_y_m"] = GUN_Y_M
        meta["t0_note"] = ("t0_pulse_ns = TOF(E_n) over (flight_path - gun_y); "
                           "hit pulse time = t0_pulse_ns + dig_t_*_ns")
        if args.t_max_us is not None:
            meta["t_max_us"] = args.t_max_us
        arrays["evt_neutron_E_eV"] = events["neutron_E_eV"].to_numpy(np.float64)
        arrays["evt_t0_pulse_ns"] = events["t0_pulse_ns"].to_numpy(np.float64)
        arrays["evt_cap_vol"] = events["capture_vol"].map(
            lambda s: CAP_VOL_CODE.get(str(s), 10)).to_numpy(np.int8)
        arrays["evt_cap_proc"] = events["capture_proc"].map(
            lambda s: CAP_PROC_CODE.get(str(s), 3)).to_numpy(np.int8)
        arrays["evt_cap_y"] = events["cap_y"].to_numpy(np.float32)

    np.savez_compressed(path, meta=json.dumps(meta), **arrays)
    size_mb = Path(path).stat().st_size / 1e6
    print(f"  -> {path}  ({len(events):,} events, "
          f"{len(digests):,} digests, {size_mb:.0f} MB)")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    files = resolve_files(args.inputs, args.max_files)
    if not files:
        print("ERROR: no input files found"); sys.exit(1)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Digesting {len(files)} files ({args.mode}) with "
          f"{args.workers} workers ...")

    work = [(f, args.mode, args.chunk_size) for f in files]
    if args.workers > 1:
        with Pool(args.workers) as pool:
            results = pool.map(digest_file, work)
    else:
        results = [digest_file(w) for w in work]
    n_err = sum(1 for r in results if "err" in r)
    if n_err:
        print(f"WARNING: {n_err}/{len(files)} files skipped")

    if args.mode == "pairs":
        for cls, name in [(0, "pool_x17.npz"), (1, "pool_ipc.npz")]:
            out = assemble_pool(results, args.mode, args, class_filter=cls)
            if out is None:
                print(f"  no events for class {cls}"); continue
            save_pool(outdir / name, *out, args.mode, args, files)
    else:
        out = assemble_pool(results, args.mode, args)
        if out is None:
            print("  no hit-bearing events found"); sys.exit(1)
        save_pool(outdir / "pool_neutron_bg.npz", *out, args.mode, args, files)

    print("Done.")


if __name__ == "__main__":
    main()
