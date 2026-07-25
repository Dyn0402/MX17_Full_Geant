#!/usr/bin/env python3
"""
analyze_trigger_thermal.py — SiPM-wall × plastics trigger optimization for the
thermal (>1 ms) window
===============================================================================
Inputs (all produced with the final 2026-07 geometry):

  --signal   pairs_thermal_trig/*.root     X17+IPC pairs, vertices from the
                                           thermal self-shielding library
  --thermal  neutrons_thermal_trig/*.root  EAR2 beam neutrons, 1 meV-2 eV
                                           (arrive in-gate by construction)
  --epi      neutrons_epi_trig/*.root      2 eV-100 keV (arrive before 1 ms;
                                           only delayed hits spill into gate)
  --mip      MIP calibration file(s) (single mu- run) — defines 1 MIP/detector

Per event and arm, the trigger observables are
  S = max single SiPM-bar edep   [MIP]   (per-channel discriminator)
  P = max plastic-bar (L/R) edep [MIP]   (per-channel discriminator)
A "leg" fires when S >= a AND P >= b; the pair tag requires >= 2 legs in
distinct arms.  Thresholds a, b are scanned; >=1-leg (singles) rates and
SiPM-only / plastic-only pair tags are recorded for comparison.

Background normalisation (EAR2 Ph3 flux, 19.5 m, TOF[ms] = 1.41/sqrt(E_n[eV])):
  thermal window [1e-3, 2] eV : 4.284e6 n/pulse   (all in-gate)
  epi window     [2, 1e5] eV  : 6.484e6 n/pulse   (in-gate iff TOF+t_hit>1 ms)

Outputs: <outdir>/trigger_scan.json + spectra/ROC figures (PDF).
"""

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import uproot

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False

# ── Constants ─────────────────────────────────────────────────────────────────
N_PULSE_THERMAL = 4.284e6    # n/pulse, [1e-3, 2] eV  (flux-file integral)
N_PULSE_EPI     = 6.484e6    # n/pulse, [2, 1e5] eV
GATE_MS         = 1.0        # gate opens 1 ms after the gamma flash
FLIGHT_TOF_MS   = 1.41       # TOF [ms] = FLIGHT_TOF_MS / sqrt(E_n [eV]), 19.5 m

SIPM_DET    = ("PlasticScint",)              # 3 mm SiPM-wall bars
PLASTIC_DET = ("BackScintL", "BackScintR")   # 25 mm plastic bars
LS_DET      = ("LiqScint_1",)

A_GRID = np.round(np.arange(0.1, 2.01, 0.1), 2)   # SiPM threshold [MIP]
B_GRID = np.round(np.arange(0.1, 2.01, 0.1), 2)   # plastic threshold [MIP]


def collect(patterns):
    files = []
    for pat in patterns:
        p = Path(pat)
        if p.is_dir():
            files += sorted(str(f) for f in p.glob("*.root"))
        else:
            files += sorted(glob.glob(pat))
    return files


def det_class_codes(det_bytes):
    """Map detType strings to codes: 0=SiPM bar, 1=plastic L, 2=plastic R,
    3=LS, -1=other.  Vectorized on the raw (possibly padded) string array."""
    det = np.asarray(det_bytes)
    code = np.full(len(det), -1, dtype=np.int8)
    def starts(prefix):
        return np.char.startswith(det.astype(str), prefix)
    code[starts("PlasticScint")] = 0
    code[starts("BackScintL")]   = 1
    code[starts("BackScintR")]   = 2
    code[starts("LiqScint_1")]   = 3
    return code


# ── MIP calibration ───────────────────────────────────────────────────────────
def mip_scale(mip_files):
    """MeV per MIP for SiPM bars / plastics / LS: Landau MPV of the per-event
    per-detector-class edep sum from a normal-incidence muon run."""
    acc = {0: {}, 1: {}, 3: {}}   # class code -> {event: sum MeV}; L+R folded to 1
    for fp in mip_files:
        with uproot.open(fp) as f:
            t = f["HitTree"].arrays(["eventID", "detType", "edep"], library="np")
        code = det_class_codes(t["detType"])
        code[code == 2] = 1
        for cls in acc:
            m = code == cls
            if not m.any():
                continue
            ev, inv = np.unique(t["eventID"][m], return_inverse=True)
            s = np.zeros(len(ev))
            np.add.at(s, inv, t["edep"][m] * 1e-6)
            for e, v in zip(ev, s):
                acc[cls][e] = acc[cls].get(e, 0.0) + v
    out = {}
    for cls, name in ((0, "sipm"), (1, "plastic"), (3, "ls")):
        v = np.array(list(acc[cls].values()))
        v = v[v > 0.05]
        if len(v) < 100:
            out[name] = None
            continue
        hi = np.quantile(v, 0.98)
        h, edges = np.histogram(v, bins=200, range=(0, hi))
        out[name] = float(0.5 * (edges[np.argmax(h)] + edges[np.argmax(h) + 1]))
    return out


# ── Event digest: per (event, arm) trigger observables, fully vectorized ──────
def _digest_chunk(t, e_n_lookup):
    """One HitTree chunk → (uev, obs) with obs (N,4,3): [S_maxbar, P_maxbar,
    LS_sum] in MeV.  e_n_lookup: None (no gate) or eventID-indexed E_n array."""
    code = det_class_codes(t["detType"])
    keep = code >= 0
    if e_n_lookup is not None:
        e_hit = e_n_lookup[t["eventID"]]
        with np.errstate(divide="ignore"):
            tof = np.where(e_hit > 0, FLIGHT_TOF_MS / np.sqrt(e_hit), np.inf)
        keep &= (tof + t["time"] * 1e-6) > GATE_MS

    code = code[keep]
    ev   = t["eventID"][keep].astype(np.int64)
    arm  = t["armID"][keep].astype(np.int64)
    edep = (t["edep"][keep] * 1e-6).astype(np.float64)     # eV → MeV
    u    = t["u"][keep]
    if len(ev) == 0:
        return np.array([], np.int64), np.zeros((0, 4, 3), np.float32)

    # channel index: SiPM bar 0..19 from u; plastics L=20 R=21; LS=22
    chan = np.full(len(ev), 22, dtype=np.int64)
    sipm = code == 0
    chan[sipm] = np.clip((u[sipm] + 250.0) // 25.0, 0, 19).astype(np.int64)
    chan[code == 1] = 20
    chan[code == 2] = 21

    # aggregate edep per (event, arm, channel)
    key = (ev * 4 + arm) * 23 + chan
    ukey, inv = np.unique(key, return_inverse=True)
    esum = np.zeros(len(ukey))
    np.add.at(esum, inv, edep)
    k_ev   = ukey // (4 * 23)
    k_arm  = (ukey // 23) % 4
    k_chan = ukey % 23

    uev, ev_inv = np.unique(k_ev, return_inverse=True)
    obs = np.zeros((len(uev), 4, 3), dtype=np.float32)
    is_s  = k_chan <= 19
    is_p  = (k_chan == 20) | (k_chan == 21)
    is_ls = k_chan == 22
    np.maximum.at(obs[:, :, 0], (ev_inv[is_s], k_arm[is_s]), esum[is_s])
    np.maximum.at(obs[:, :, 1], (ev_inv[is_p], k_arm[is_p]), esum[is_p])
    np.add.at(obs[:, :, 2], (ev_inv[is_ls], k_arm[is_ls]), esum[is_ls])
    return uev, obs


def digest_file(fp, gate=False, truth=False):
    """Returns (obs, aux, n_events_total): obs (N,4,3) float32 per surviving
    event & arm [S_maxbar, P_maxbar, LS_sum] in MeV; aux: per-row eventID and
    (if truth) event_type / min pair KE.  Streams HitTree in chunks."""
    ev_branches = ["eventID", "event_type", "neutron_E_eV"]
    if truth:
        ev_branches += ["em_ke", "ep_ke"]
    parts = []
    with uproot.open(fp) as f:
        nev = f["EventTree"].num_entries
        et = f["EventTree"].arrays(ev_branches, library="np")
        e_n_lookup = None
        if gate and nev:
            e_n_lookup = np.zeros(int(et["eventID"].max()) + 1)
            e_n_lookup[et["eventID"]] = et["neutron_E_eV"]
        if f["HitTree"].num_entries:
            for t in f["HitTree"].iterate(
                    ["eventID", "armID", "detType", "edep", "time", "u"],
                    library="np", step_size="200 MB"):
                parts.append(_digest_chunk(t, e_n_lookup))

    parts = [(u_, o_) for (u_, o_) in parts if len(u_)]
    if not parts:
        return (np.zeros((0, 4, 3), np.float32),
                {"eventID": np.array([], np.int64)}, nev)
    # merge chunks (an event may straddle a boundary): max for S/P, add for LS
    all_ev  = np.concatenate([u_ for u_, _ in parts])
    all_obs = np.concatenate([o_ for _, o_ in parts])
    uev, inv = np.unique(all_ev, return_inverse=True)
    obs = np.zeros((len(uev), 4, 3), dtype=np.float32)
    for a in range(4):
        np.maximum.at(obs[:, a, 0], inv, all_obs[:, a, 0])
        np.maximum.at(obs[:, a, 1], inv, all_obs[:, a, 1])
        np.add.at(obs[:, a, 2], inv, all_obs[:, a, 2])

    aux = {"eventID": uev}
    if truth and nev:
        idx = np.zeros(int(et["eventID"].max()) + 1, dtype=np.int64)
        idx[et["eventID"]] = np.arange(nev)
        aux["event_type"] = et["event_type"][idx[uev]]
        aux["ke_min"] = np.minimum(et["em_ke"], et["ep_ke"])[idx[uev]]
    return obs, aux, nev


# ── Trigger scan (vectorized over events, looped over the 20×20 grid) ─────────
def scan(obs, mip_sipm, mip_plastic):
    """obs: (N,4,3) MeV.  Returns leg1/leg2 counts on the (a,b) grid and
    SiPM-only / plastic-only >=2-arm counts on the 1-D grids."""
    S = obs[:, :, 0] / mip_sipm      # (N,4) [MIP]
    P = obs[:, :, 1] / mip_plastic
    nA, nB = len(A_GRID), len(B_GRID)
    n1 = np.zeros((nA, nB), np.int64)
    n2 = np.zeros((nA, nB), np.int64)
    n2m = np.zeros((nA, nB), np.int64)   # 2 SiPM arms + >=1 plastic confirm
    for ia, a in enumerate(A_GRID):
        sa = S >= a
        two_sipm = sa.sum(axis=1) >= 2
        for ib, b in enumerate(B_GRID):
            legs = sa & (P >= b)
            nlegs = legs.sum(axis=1)
            n1[ia, ib] = (nlegs >= 1).sum()
            n2[ia, ib] = (nlegs >= 2).sum()
            n2m[ia, ib] = (two_sipm & (nlegs >= 1)).sum()
    n2s = np.array([((S >= a).sum(axis=1) >= 2).sum() for a in A_GRID])
    n2p = np.array([((P >= b).sum(axis=1) >= 2).sum() for b in B_GRID])
    return n1, n2, n2m, n2s, n2p


def _sig_digest(fp):
    obs, aux, nev = digest_file(fp, truth=True)
    with uproot.open(fp) as f:
        ty = f["EventTree"]["event_type"].array(library="np")
    return (obs, aux.get("event_type", np.array([], np.int32)),
            aux.get("ke_min", np.array([])),
            int((ty == 0).sum()), int((ty == 1).sum()))


def _bg_digest(fp_gate):
    fp, gate = fp_gate
    obs, _, nev = digest_file(fp, gate=gate)
    return obs, nev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signal",  nargs="+", required=True)
    ap.add_argument("--thermal", nargs="+", required=True)
    ap.add_argument("--epi",     nargs="+", default=[])
    ap.add_argument("--mip",     nargs="+", default=[])
    ap.add_argument("--mip-mev", default=None,
                    help="explicit 'sipm,plastic[,ls]' MeV per MIP; skips --mip. "
                         "Geometry-invariant to the capsule flip, so the "
                         "valve-first mu- values stay valid nose-first.")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("-o", "--outdir", default="analysis/trigger_thermal")
    ap.add_argument("--max-files", type=int, default=None)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sig_files = collect(args.signal)[: args.max_files]
    th_files  = collect(args.thermal)[: args.max_files]
    ep_files  = collect(args.epi)[: args.max_files]
    for name, fl in (("signal", sig_files), ("thermal", th_files)):
        if not fl:
            print(f"ERROR: no {name} files", file=sys.stderr)
            sys.exit(1)

    if args.mip_mev:
        v = [float(x) for x in args.mip_mev.split(",")]
        mips = {"sipm": v[0], "plastic": v[1], "ls": v[2] if len(v) > 2 else None}
        print("MIP calibration (explicit --mip-mev)")
    else:
        mip_files = collect(args.mip)
        if not mip_files:
            print("ERROR: no --mip files and no --mip-mev given", file=sys.stderr)
            sys.exit(1)
        print(f"MIP calibration from {len(mip_files)} file(s) ...")
        mips = mip_scale(mip_files)
        if not mips["sipm"] or not mips["plastic"]:
            print("ERROR: MIP MPV extraction failed", file=sys.stderr)
            sys.exit(1)
    print(f"  1 MIP = SiPM bar {mips['sipm']*1e3:.0f} keV | "
          f"plastic {mips['plastic']:.2f} MeV | "
          f"LS {mips['ls'] if mips['ls'] else 'n/a'} MeV")

    from concurrent.futures import ProcessPoolExecutor

    def pmap(fn, items):
        if args.workers > 1 and len(items) > 1:
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                return list(pool.map(fn, items))
        return [fn(x) for x in items]

    # ── Signal ────────────────────────────────────────────────────────────
    sr = pmap(_sig_digest, sig_files)
    sig_obs   = np.concatenate([r[0] for r in sr])
    sig_type  = np.concatenate([r[1] for r in sr])
    sig_kemin = np.concatenate([r[2] for r in sr])
    n_gen = {0: sum(r[3] for r in sr), 1: sum(r[4] for r in sr)}
    print(f"Signal: {n_gen[0]:,} X17 + {n_gen[1]:,} IPC generated; "
          f"{len(sig_obs):,} events with hits")

    # ── Backgrounds ───────────────────────────────────────────────────────
    def run_bg(files, gate, label):
        if not files:
            return np.zeros((0, 4, 3), np.float32), 0
        rr = pmap(_bg_digest, [(fp, gate) for fp in files])
        obs = np.concatenate([r[0] for r in rr])
        n_tot = sum(r[1] for r in rr)
        print(f"{label}: {n_tot:,} neutrons, {len(obs):,} events with "
              f"in-gate detector hits")
        return obs, n_tot

    th_obs, n_th = run_bg(th_files, gate=False, label="thermal bg")
    ep_obs, n_ep = (run_bg(ep_files, gate=True, label="epi spill-in")
                    if ep_files else (np.zeros((0, 4, 3), np.float32), 0))

    # ── Scans ─────────────────────────────────────────────────────────────
    res = {"mip_MeV": mips, "a_grid": A_GRID.tolist(), "b_grid": B_GRID.tolist(),
           "n_pulse_thermal": N_PULSE_THERMAL, "n_pulse_epi": N_PULSE_EPI,
           "gate_ms": GATE_MS}

    classes = [
        ("x17",     sig_obs[sig_type == 0], n_gen[0]),
        ("ipc",     sig_obs[sig_type == 1], n_gen[1]),
        ("ipc_hi",  sig_obs[(sig_type == 1) & (sig_kemin > 5.0)],
                    int(((sig_type == 1) & (sig_kemin > 5.0)).sum())),
        ("thermal", th_obs, n_th),
        ("epi",     ep_obs, n_ep),
    ]
    for label, obs, norm in classes:
        if norm == 0:
            continue
        n1, n2, n2m, n2s, n2p = scan(obs, mips["sipm"], mips["plastic"])
        res[label] = {"n": int(norm),
                      "leg1": (n1 / norm).tolist(),
                      "leg2": (n2 / norm).tolist(),
                      "leg2_confirm1": (n2m / norm).tolist(),
                      "leg2_sipm_only": (n2s / norm).tolist(),
                      "leg2_plastic_only": (n2p / norm).tolist()}

    with open(outdir / "trigger_scan.json", "w") as f:
        json.dump(res, f, indent=1)
    print(f"\nWrote {outdir}/trigger_scan.json")

    # ── Headline table ────────────────────────────────────────────────────
    e2  = np.array(res["x17"]["leg2"])
    ei2 = np.array(res["ipc"]["leg2"])
    bg1 = np.array(res["thermal"]["leg1"]) * N_PULSE_THERMAL
    bg2 = np.array(res["thermal"]["leg2"]) * N_PULSE_THERMAL
    if "epi" in res:
        bg1 = bg1 + np.array(res["epi"]["leg1"]) * N_PULSE_EPI
        bg2 = bg2 + np.array(res["epi"]["leg2"]) * N_PULSE_EPI
    em2  = np.array(res["x17"]["leg2_confirm1"])
    bgm2 = np.array(res["thermal"]["leg2_confirm1"]) * N_PULSE_THERMAL
    if "epi" in res:
        bgm2 = bgm2 + np.array(res["epi"]["leg2_confirm1"]) * N_PULSE_EPI
    print("\n a[MIP] b[MIP]  eff2(X17) eff2c(X17) eff2(IPC)  "
          "pairBG/pulse  confBG/pulse  legBG/pulse")
    show = (0.2, 0.5, 1.0)
    for ia, a in enumerate(A_GRID):
        for ib, b in enumerate(B_GRID):
            if a in show and b in show:
                print(f"  {a:4.1f}  {b:4.1f}   {e2[ia,ib]:.3f}     "
                      f"{em2[ia,ib]:.3f}     {ei2[ia,ib]:.3f}     "
                      f"{bg2[ia,ib]:.3e}    {bgm2[ia,ib]:.3e}    "
                      f"{bg1[ia,ib]:.3e}")

    # ── Figures ───────────────────────────────────────────────────────────
    if not HAVE_MPL:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, idx, mip, name in ((axes[0], 0, mips["sipm"], "SiPM bar max"),
                               (axes[1], 1, mips["plastic"], "plastic bar max")):
        bins = np.linspace(0, 6, 121)
        for obs, lab, sty in ((sig_obs[sig_type == 0], "X17", {}),
                              (sig_obs[sig_type == 1], "IPC", {"ls": "--"}),
                              (th_obs, "thermal bg", {"ls": ":"})):
            v = (obs[:, :, idx] / mip).ravel()
            v = v[v > 0]
            if len(v):
                ax.hist(v, bins=bins, histtype="step", density=True,
                        label=lab, **sty)
        ax.set_xlabel(f"{name} [MIP]")
        ax.set_yscale("log")
        ax.legend()
    fig.suptitle("Per-arm trigger observables, thermal (>1 ms) window")
    fig.tight_layout()
    fig.savefig(outdir / "trigger_spectra.pdf")

    fig, ax = plt.subplots(figsize=(6, 4.5))
    for ia, a in enumerate(A_GRID):
        if a in (0.2, 0.3, 0.5, 1.0):
            ax.plot(bg2[ia], e2[ia], marker="o", ms=3,
                    label=f"2 full legs, SiPM ≥ {a:.1f} MIP")
            ax.plot(bgm2[ia], em2[ia], marker="s", ms=3, ls="--",
                    label=f"2 SiPM + 1 confirm, SiPM ≥ {a:.1f} MIP")
    ax.set_xscale("log")
    ax.set_xlabel("pair-tag background / pulse (>1 ms gate)")
    ax.set_ylabel("X17 pair-tag efficiency")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "trigger_roc.pdf")
    print(f"Figures → {outdir}/trigger_spectra.pdf, trigger_roc.pdf")


if __name__ == "__main__":
    main()
