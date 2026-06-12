#!/usr/bin/env python3
"""
analyze_mev_captures.py
Full-range capture scan of run B-full (neutron-beam mode, 1 meV - 100 MeV).

The MeV analogue of analyze_thermal_captures.py, with one crucial difference:
the radiative rate is taken from the DIRECT 3He(n,g) counts.  At full-range
statistics (5e8 events) the direct channel has hundreds of events, and the
thermal sigma_ng/sigma_np = 1.013e-8 used in the sub-keV analysis is wrong by
3-4 orders of magnitude at MeV energies (the ratio rises to ~1e-4; see
ENDF MT102/MT103 in data/He3.h5).  The (n,p) x thermal-ratio product is kept
only as the sub-keV opaque-gas floor.

Per-pulse normalisation: the generator samples flux_n_pulse_NOisolet_100bpd
(data/fluxEAR2-Ph3_in_different_units.root) over 1 meV - 100 MeV; that
histogram integrates to 2.263e7 n/pulse (its < 1 keV part reproduces the
7.31e6 sub-keV anchor).  Do NOT reuse the sub-keV anchor as the total.

Usage (lxplus):
    source scripts/setup_lxplus.sh
    python3 scripts/analyze_mev_captures.py \
        /eos/experiment/ntof/data/x17/full_sim/neutrons_fullrange/ \
        -o mev_captures --workers 8
Figures + derived rate tables: scripts/make_mev_report_figures.py (reads the
.npz written here).
"""

import argparse
import glob
import json
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import uproot

# ── Constants ─────────────────────────────────────────────────────────────────
N_PER_PULSE_DEFAULT = 2.2628e7    # flux-file integral, 1 meV - 100 MeV
PULSES_PER_DAY      = 1.929e4     # 7e12 ppp (results_3He header)
SIGMA_RATIO_THERMAL = 54e-6 / 5333.0   # 1.013e-8; sub-keV floor ONLY
ALPHA_IPC_DEFAULT   = 2.1e-3      # IPC per radiative capture (Alberto)
BR_X17_DEFAULT      = 0.025       # X17 per IPC pair (NOT per capture)

# log10(E/eV) binning: 1 meV - 100 MeV, 10 bins/decade
LOGE_MIN, LOGE_MAX, NBINS = -3.0, 8.0, 110
LOGE_EDGES = np.linspace(LOGE_MIN, LOGE_MAX, NBINS + 1)

# Alberto thin-target table (results_3He), per decade:
# decade lower edge [eV] -> He3-captures per incident neutron
# (= He3-captures/pulse / neutrons/pulse, so his flux assumption divides out)
ALBERTO_THIN_TARGET = {
    1e-2: 4.37e0/3.11e6, 1e-1: 9.60e-1/1.68e6, 1e0: 2.33e-1/1.39e6,
    1e1: 1.05e-1/1.57e6, 1e2: 1.05e-1/1.78e6,  1e3: 2.90e-1/2.09e6,
    1e4: 1.32e0/3.03e6,  1e5: 1.50e1/9.89e6,   1e6: 1.77e1/6.75e6,
    1e7: 1.56e0/1.64e6,
}

SCINT_VOLS = ("PlasticScint", "LiqScint_1", "LiqScint_2",
              "BackScintL", "BackScintR")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="*",
                   help="Run-B-full ROOT files, glob(s), or a directory")
    p.add_argument("-o", "--output", default="mev_captures",
                   help="Output prefix (.npz, .json)")
    p.add_argument("--n-per-pulse", type=float, default=N_PER_PULSE_DEFAULT,
                   help="Beam neutrons/pulse over the sampled energy range "
                        "(default: flux-file integral 1 meV-100 MeV)")
    p.add_argument("--alpha-ipc",   type=float, default=ALPHA_IPC_DEFAULT)
    p.add_argument("--br-x17",      type=float, default=BR_X17_DEFAULT)
    p.add_argument("--workers",     type=int,   default=4)
    p.add_argument("--max-files",   type=int,   default=None)
    p.add_argument("--step-size",   type=int,   default=2_000_000,
                   help="EventTree rows per chunk")
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
    if max_files:
        files = files[:max_files]
    return files


# ── Per-file scan (runs in worker process) ────────────────────────────────────
def scan_file(args):
    try:
        return {"ok": True, "acc": _scan_file_impl(args)}
    except Exception as e:
        print(f"[SKIP] {Path(args[0]).name}: {e.__class__.__name__}: "
              f"{str(e)[:80]}", flush=True)
        return {"ok": False, "path": args[0], "err": str(e)[:200]}


def _scan_file_impl(args):
    path, step_size = args
    acc = {
        "n_events": 0,
        "budget": Counter(),                       # (vol, proc) -> count
        "h_all":      np.zeros(NBINS, dtype=np.int64),
        "h_gas_np":   np.zeros(NBINS, dtype=np.int64),
        "h_gas_rad":  np.zeros(NBINS, dtype=np.int64),
        "h_wall_al":  np.zeros(NBINS, dtype=np.int64),
        "h_scint_h":  np.zeros(NBINS, dtype=np.int64),
        "h_cfrp":     np.zeros(NBINS, dtype=np.int64),
        "h_plastic":  np.zeros(NBINS, dtype=np.int64),
        "rad_events": [],                          # E_n of each direct (n,g)
    }
    branches = ["neutron_E_eV", "capture_vol", "capture_proc"]
    with uproot.open(path) as f:
        tree = f["EventTree"]
        for chunk in tree.iterate(branches, step_size=step_size, library="np"):
            E    = chunk["neutron_E_eV"]
            # Materialise awkward string arrays -> numpy object arrays of str
            # (see docs/uproot_awkward_concat_pitfall.md)
            vol  = np.asarray(chunk["capture_vol"],  dtype=object)
            proc = np.asarray(chunk["capture_proc"], dtype=object)
            acc["n_events"] += len(E)

            u_v, inv_v = np.unique(vol,  return_inverse=True)
            u_p, inv_p = np.unique(proc, return_inverse=True)
            joint = inv_v.astype(np.int64) * len(u_p) + inv_p
            cnt = np.bincount(joint, minlength=len(u_v) * len(u_p))
            for code in np.nonzero(cnt)[0]:
                v, pr = u_v[code // len(u_p)], u_p[code % len(u_p)]
                acc["budget"][(v or "(escaped)", pr or "(none)")] += int(cnt[code])

            logE = np.log10(np.clip(E, 10**LOGE_MIN, 10**LOGE_MAX))
            acc["h_all"] += np.histogram(logE, LOGE_EDGES)[0]

            in_gas  = vol == "He3Gas"
            gas_np  = in_gas & (proc == "neutronInelastic")   # 3He(n,p)t
            gas_rad = in_gas & (proc == "nCapture")           # 3He(n,g): signal
            wall_al = (vol == "He3Cap_Al") & (proc == "nCapture")
            scint_h = np.isin(vol, SCINT_VOLS) & (proc == "nCapture")
            cfrp    = (vol == "He3Cap_CFRP") & (proc == "nCapture")
            plastic = (vol == "PlasticScint") & (proc == "nCapture")

            acc["h_gas_np"]  += np.histogram(logE[gas_np],  LOGE_EDGES)[0]
            acc["h_gas_rad"] += np.histogram(logE[gas_rad], LOGE_EDGES)[0]
            acc["h_wall_al"] += np.histogram(logE[wall_al], LOGE_EDGES)[0]
            acc["h_scint_h"] += np.histogram(logE[scint_h], LOGE_EDGES)[0]
            acc["h_cfrp"]    += np.histogram(logE[cfrp],    LOGE_EDGES)[0]
            acc["h_plastic"] += np.histogram(logE[plastic], LOGE_EDGES)[0]
            acc["rad_events"].extend(float(e) for e in E[gas_rad])

    print(f"[done] {Path(path).name}  ({acc['n_events']:,} events, "
          f"{len(acc['rad_events'])} radiative)", flush=True)
    return acc


def merge_accs(accs):
    out = accs[0]
    for a in accs[1:]:
        out["n_events"] += a["n_events"]
        out["budget"]   += a["budget"]
        for k in ("h_all", "h_gas_np", "h_gas_rad", "h_wall_al", "h_scint_h",
                  "h_cfrp", "h_plastic"):
            out[k] += a[k]
        out["rad_events"] += a["rad_events"]
    return out


# ── Summary + persistence ─────────────────────────────────────────────────────
def decade_table(acc, w, alpha_ipc, br_x17):
    """Per-decade rows.  Radiative rate = DIRECT counts (Poisson stats);
    (n,p) x thermal ratio is reported only as the sub-keV floor."""
    rows = []
    for d in range(int(LOGE_MAX - LOGE_MIN)):
        lo, hi = d * 10, (d + 1) * 10
        elo = 10.0 ** (LOGE_MIN + d)
        n_beam = int(acc["h_all"][lo:hi].sum())
        n_np   = int(acc["h_gas_np"][lo:hi].sum())
        n_rad  = int(acc["h_gas_rad"][lo:hi].sum())
        rad_pp = n_rad * w
        rows.append({
            "E_lo_eV": elo, "E_hi_eV": elo * 10,
            "n_beam_sim": n_beam, "n_np_sim": n_np, "n_rad_direct": n_rad,
            "beam_per_pulse": n_beam * w,
            "np_frac_of_beam": n_np / n_beam if n_beam else 0.0,
            "rad_per_pulse": rad_pp,
            "rad_per_n": n_rad / n_beam if n_beam else 0.0,
            "subkev_floor_per_pulse": n_np * SIGMA_RATIO_THERMAL * w,
            "ipc_per_pulse": rad_pp * alpha_ipc,
            "x17_per_day": rad_pp * alpha_ipc * br_x17 * PULSES_PER_DAY,
            "thin_target_per_n": ALBERTO_THIN_TARGET.get(elo),
        })
    return rows


def print_summary(acc, args):
    w = args.n_per_pulse / acc["n_events"]
    total = acc["n_events"]
    print(f"\n{'='*80}")
    print(f"  Scanned events : {total:,}")
    print(f"  Weight/neutron : {w:.3e} per pulse "
          f"({args.n_per_pulse:.4g} n/pulse / {total:.3g} simulated)")
    print(f"  Pulses/day     : {PULSES_PER_DAY:.4g}")
    print(f"{'='*80}")

    print("\n--- Terminal-interaction budget (top 20) ---")
    for (volp, count) in acc["budget"].most_common(20):
        print(f"  {volp[0]:<24s} {volp[1]:<20s} {count:>12,}  "
              f"({count/total:.3e}/n)")

    n_rad = len(acc["rad_events"])
    n_np  = int(acc["h_gas_np"].sum())
    rad_pp = n_rad * w
    print(f"\n--- 3He radiative capture (X17 channel, DIRECT counts) ---")
    print(f"  3He(n,p)t          : {n_np:,}  ({n_np/total:.4f}/n)")
    print(f"  3He(n,g)4He direct : {n_rad}  ({n_rad/total:.3e}/n; "
          f"stat err {1/np.sqrt(max(n_rad,1)):.1%})")
    print(f"  rad capt / pulse   : {rad_pp:.2f}   / day: {rad_pp*PULSES_PER_DAY:.3e}")
    print(f"  -> IPC / pulse     : {rad_pp*args.alpha_ipc:.3e}   "
          f"/ day: {rad_pp*args.alpha_ipc*PULSES_PER_DAY:.1f}")
    print(f"  -> X17 / day       : {rad_pp*args.alpha_ipc*args.br_x17*PULSES_PER_DAY:.1f}"
          f"   (BR={args.br_x17} per IPC pair)")

    print("\n--- Per-decade (direct counts) ---")
    print(f"  {'E range':<18s}{'beam/pulse':>12s}{'(n,p)/beam':>11s}"
          f"{'N_rad':>7s}{'capt/pulse':>12s}{'X17/day':>9s}{'G4/table':>9s}")
    for r in decade_table(acc, w, args.alpha_ipc, args.br_x17):
        tt = r["thin_target_per_n"]
        g4tab = f"{r['rad_per_n']/tt:.2f}" if (tt and r["n_rad_direct"]) else "-"
        print(f"  {r['E_lo_eV']:.0e}-{r['E_hi_eV']:.0e} eV"
              f"{r['beam_per_pulse']:>12.3e}{r['np_frac_of_beam']:>11.3f}"
              f"{r['n_rad_direct']:>7d}{r['rad_per_pulse']:>12.3f}"
              f"{r['x17_per_day']:>9.3f}{g4tab:>9s}")


def save_results(acc, args):
    rad = np.array(sorted(acc["rad_events"]), dtype=float).reshape(-1, 1)
    np.savez(args.output + ".npz",
             n_events=acc["n_events"],
             loge_edges=LOGE_EDGES,
             h_all=acc["h_all"], h_gas_np=acc["h_gas_np"],
             h_gas_rad=acc["h_gas_rad"], h_wall_al=acc["h_wall_al"],
             h_scint_h=acc["h_scint_h"], h_cfrp=acc["h_cfrp"],
             h_plastic=acc["h_plastic"],
             rad_events=rad,
             budget_keys=np.array([f"{v}|{p}" for (v, p) in acc["budget"]]),
             budget_counts=np.array(list(acc["budget"].values())),
             n_per_pulse=args.n_per_pulse,
             alpha_ipc=args.alpha_ipc, br_x17=args.br_x17,
             pulses_per_day=PULSES_PER_DAY)
    w = args.n_per_pulse / acc["n_events"]
    n_rad = len(acc["rad_events"])
    summary = {
        "n_events": acc["n_events"], "n_per_pulse": args.n_per_pulse,
        "pulses_per_day": PULSES_PER_DAY, "weight_per_neutron": w,
        "n_gas_np": int(acc["h_gas_np"].sum()),
        "n_gas_rad_direct": n_rad,
        "rad_per_pulse": n_rad * w,
        "ipc_per_pulse": n_rad * w * args.alpha_ipc,
        "x17_per_day": n_rad * w * args.alpha_ipc * args.br_x17 * PULSES_PER_DAY,
        "decades": decade_table(acc, w, args.alpha_ipc, args.br_x17),
        "budget": {f"{v}|{p}": int(c) for (v, p), c in
                   acc["budget"].most_common()},
        "skipped_files": acc.get("skipped_files", []),
    }
    with open(args.output + ".json", "w") as jf:
        json.dump(summary, jf, indent=1)
    print(f"\nSaved: {args.output}.npz  {args.output}.json")
    print("Figures + report tables: python3 scripts/make_mev_report_figures.py")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    files = resolve_files(args.inputs, args.max_files)
    if not files:
        print("ERROR: no input files found"); sys.exit(1)
    print(f"Scanning {len(files)} files with {args.workers} workers ...")

    work = [(f, args.step_size) for f in files]
    if args.workers > 1:
        with Pool(args.workers) as pool:
            results = pool.map(scan_file, work)
    else:
        results = [scan_file(wk) for wk in work]
    accs    = [r["acc"] for r in results if r["ok"]]
    skipped = [(r["path"], r["err"]) for r in results if not r["ok"]]
    if skipped:
        print(f"\nWARNING: skipped {len(skipped)}/{len(files)} unreadable files:")
        for p, _ in skipped:
            print(f"  {Path(p).name}")
    if not accs:
        print("ERROR: no readable files"); sys.exit(1)
    acc = merge_accs(accs)
    acc["skipped_files"] = [Path(p).name for p, _ in skipped]

    print_summary(acc, args)
    save_results(acc, args)


if __name__ == "__main__":
    main()
