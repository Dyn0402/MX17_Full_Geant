#!/usr/bin/env python3
"""
analyze_thermal_captures.py
Thermal-statistics analysis of run B (neutron-beam mode, event_type 2).

Question: how many radiative captures 3He(n,g)4He -- i.e. X17/IPC production
events -- does the EAR2 beam give per pulse as a function of neutron energy,
with full transport (self-shielding, scattering) included?

Reads the EventTree of every job file and accumulates:
  - terminal-interaction budget per (volume, process)
  - log10(E_n) histograms: all primaries, 3He(n,p)t, 3He(n,g) [radiative],
    He3Cap_Al nCapture (wall), H-capture in scintillator volumes
  - every individual radiative-capture event (E_n, cap_y) -- expected O(10)
  - absorption depth (cap_y) of (n,p)t events in three E_n slices

Per-pulse normalisation: each simulated neutron represents
(n_per_pulse / N_simulated) real neutrons; the generator already sampled the
EAR2 flux shape and the energy-dependent transverse profile, so the weighted
histograms ARE per-pulse rates.

Usage (lxplus):
    source scripts/setup_lxplus.sh
    python3 scripts/analyze_thermal_captures.py \
        /eos/experiment/ntof/data/x17/full_sim/neutrons_subkev/ \
        -o thermal_captures --workers 8
    # plots from a previous scan:
    python3 scripts/analyze_thermal_captures.py --plot-only thermal_captures.npz
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
N_PER_PULSE_DEFAULT = 7.31e6      # EAR2 neutrons/pulse below 1 keV (flux file)
SIGMA_RATIO         = 54e-6 / 5333.0   # s_ng / s_np at thermal = 1.013e-8
ALPHA_IPC_DEFAULT   = 2.1e-3      # IPC per radiative capture (Alberto's table)
BR_X17_DEFAULT      = 0.025       # assumed X17 per radiative capture

# log10(E/eV) binning: 1 meV - 1 keV, 10 bins/decade
LOGE_MIN, LOGE_MAX, NBINS = -3.0, 3.0, 60
LOGE_EDGES = np.linspace(LOGE_MIN, LOGE_MAX, NBINS + 1)

# Alberto's thin-target radiative captures per incident neutron, per decade
# (calculation_tables/results_3He; see docs/he3_self_shielding_note.md)
ALBERTO_THIN_TARGET = {  # (decade lower edge in eV) -> rad capt / n
    1e-2: 1.41e-6, 1e-1: 5.71e-7, 1e0: 1.68e-7, 1e1: 6.69e-8, 1e2: 5.90e-8,
}
ANALYTIC_CAP = {  # optically-thick cap s_ng/s_np per decade
    1e-2: 1.01e-8, 1e-1: 1.01e-8, 1e0: 1.01e-8, 1e1: 1.00e-8, 1e2: 7.9e-9,
}

SCINT_VOLS = ("PlasticScint", "LiqScint_1", "LiqScint_2",
              "BackScintL", "BackScintR")

# Absorption-depth slices (eV) for the cap_y validation histogram
DEPTH_SLICES = [(1e-3, 1e-1), (1e-1, 1e1), (1e1, 1e3)]
CAPY_EDGES = np.linspace(-30.0, 30.0, 121)  # mm, beam axis through gas


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="*",
                   help="Run-B ROOT files, glob(s), or a directory")
    p.add_argument("-o", "--output", default="thermal_captures",
                   help="Output prefix (.npz, .json, .pdf)")
    p.add_argument("--n-per-pulse", type=float, default=N_PER_PULSE_DEFAULT)
    p.add_argument("--alpha-ipc",   type=float, default=ALPHA_IPC_DEFAULT)
    p.add_argument("--br-x17",      type=float, default=BR_X17_DEFAULT)
    p.add_argument("--workers",     type=int,   default=4)
    p.add_argument("--max-files",   type=int,   default=None)
    p.add_argument("--step-size",   type=int,   default=2_000_000,
                   help="EventTree rows per chunk")
    p.add_argument("--plot", action="store_true",
                   help="Also produce the PDF figure after scanning")
    p.add_argument("--plot-only", metavar="NPZ", default=None,
                   help="Skip the scan; plot from an existing .npz")
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
    path, step_size = args
    acc = {
        "n_events": 0,
        "budget": Counter(),                       # (vol, proc) -> count
        "h_all":      np.zeros(NBINS, dtype=np.int64),
        "h_gas_np":   np.zeros(NBINS, dtype=np.int64),
        "h_gas_rad":  np.zeros(NBINS, dtype=np.int64),
        "h_wall_al":  np.zeros(NBINS, dtype=np.int64),
        "h_scint_h":  np.zeros(NBINS, dtype=np.int64),
        "rad_events": [],                          # (E_eV, cap_y) per event
        "h_depth": [np.zeros(len(CAPY_EDGES) - 1, dtype=np.int64)
                    for _ in DEPTH_SLICES],
    }
    branches = ["neutron_E_eV", "capture_vol", "capture_proc", "cap_y"]
    with uproot.open(path) as f:
        tree = f["EventTree"]
        for chunk in tree.iterate(branches, step_size=step_size, library="np"):
            E    = chunk["neutron_E_eV"]
            # Materialise awkward string arrays -> numpy object arrays of str
            # (see docs/uproot_awkward_concat_pitfall.md)
            vol  = np.asarray(chunk["capture_vol"],  dtype=object)
            proc = np.asarray(chunk["capture_proc"], dtype=object)
            capy = chunk["cap_y"]
            n = len(E)
            acc["n_events"] += n

            # (vol, proc) budget via joint integer codes
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
            gas_np  = in_gas & (proc == "neutronInelastic")
            gas_rad = in_gas & (proc == "nCapture")
            wall_al = (vol == "He3Cap_Al") & (proc == "nCapture")
            scint_h = np.isin(vol, SCINT_VOLS) & (proc == "nCapture")

            acc["h_gas_np"]  += np.histogram(logE[gas_np],  LOGE_EDGES)[0]
            acc["h_gas_rad"] += np.histogram(logE[gas_rad], LOGE_EDGES)[0]
            acc["h_wall_al"] += np.histogram(logE[wall_al], LOGE_EDGES)[0]
            acc["h_scint_h"] += np.histogram(logE[scint_h], LOGE_EDGES)[0]

            for E_rad, y_rad in zip(E[gas_rad], capy[gas_rad]):
                acc["rad_events"].append((float(E_rad), float(y_rad)))

            for i, (elo, ehi) in enumerate(DEPTH_SLICES):
                m = gas_np & (E >= elo) & (E < ehi)
                acc["h_depth"][i] += np.histogram(capy[m], CAPY_EDGES)[0]
    return acc


def merge_accs(accs):
    out = accs[0]
    for a in accs[1:]:
        out["n_events"] += a["n_events"]
        out["budget"]   += a["budget"]
        for k in ("h_all", "h_gas_np", "h_gas_rad", "h_wall_al", "h_scint_h"):
            out[k] += a[k]
        out["rad_events"] += a["rad_events"]
        for i in range(len(DEPTH_SLICES)):
            out["h_depth"][i] += a["h_depth"][i]
    return out


# ── Summary + persistence ─────────────────────────────────────────────────────
def decade_table(acc, w):
    """Aggregate the 10-bins/decade histograms into per-decade rows."""
    rows = []
    for d in range(int(LOGE_MAX - LOGE_MIN)):
        lo, hi = d * 10, (d + 1) * 10
        elo = 10.0 ** (LOGE_MIN + d)
        n_beam = int(acc["h_all"][lo:hi].sum())
        n_np   = int(acc["h_gas_np"][lo:hi].sum())
        n_rad  = int(acc["h_gas_rad"][lo:hi].sum())
        rows.append({
            "E_lo_eV": elo, "E_hi_eV": elo * 10,
            "n_beam_sim": n_beam, "n_np_sim": n_np, "n_rad_sim": n_rad,
            "beam_per_pulse": n_beam * w,
            "np_per_pulse":   n_np * w,
            "rad_per_pulse_direct": n_rad * w,
            "rad_per_pulse_npxbr":  n_np * SIGMA_RATIO * w,
            "rad_per_n_direct": n_rad / n_beam if n_beam else 0.0,
            "rad_per_n_npxbr":  n_np * SIGMA_RATIO / n_beam if n_beam else 0.0,
            "thin_target_per_n": ALBERTO_THIN_TARGET.get(elo),
            "cap_per_n": ANALYTIC_CAP.get(elo),
        })
    return rows


def print_summary(acc, args):
    w = args.n_per_pulse / acc["n_events"]
    print(f"\n{'='*78}")
    print(f"  Scanned events : {acc['n_events']:,}")
    print(f"  Weight/neutron : {w:.3e} per pulse "
          f"({args.n_per_pulse:.3g} n/pulse / {acc['n_events']:.3g} simulated)")
    print(f"{'='*78}")

    print("\n--- Terminal-interaction budget (top 20) ---")
    total = acc["n_events"]
    for (volp, count) in acc["budget"].most_common(20):
        print(f"  {volp[0]:<24s} {volp[1]:<20s} {count:>12,}  "
              f"({count/total:.3e}/n)")

    n_rad = len(acc["rad_events"])
    n_np  = int(acc["h_gas_np"].sum())
    print(f"\n--- 3He radiative capture (X17 channel) ---")
    print(f"  3He(n,p)t          : {n_np:,}  ({n_np/total:.4f}/n)")
    print(f"  3He(n,g)4He direct : {n_rad}  ({n_rad/total:.3e}/n)")
    if n_np:
        ratio = n_rad / n_np
        print(f"  ratio ng/np        : {ratio:.3e}  "
              f"(expected s_ng/s_np = {SIGMA_RATIO:.3e})")
    print(f"  rad capt / pulse   : direct {n_rad*w:.3e}   "
          f"(n,p)xBR {n_np*SIGMA_RATIO*w:.3e}")
    rad_pp = n_np * SIGMA_RATIO * w
    print(f"  -> IPC / pulse     : {rad_pp*args.alpha_ipc:.3e}  "
          f"(alpha_IPC={args.alpha_ipc})")
    print(f"  -> X17 / pulse     : {rad_pp*args.br_x17:.3e}  "
          f"(BR={args.br_x17})")
    print(f"  [thin-target table : 1.21e-2 IPC/pulse; "
          f"self-shielding note: 1.9e-4]")

    print("\n--- Per-decade ---")
    hdr = (f"  {'E range':<18s}{'beam/pulse':>12s}{'(n,p)/n':>10s}"
           f"{'rad/n direct':>14s}{'rad/n (np)xBR':>14s}{'thin-tgt':>10s}{'cap':>10s}")
    print(hdr)
    for r in decade_table(acc, w):
        tt = f"{r['thin_target_per_n']:.2e}" if r['thin_target_per_n'] else "-"
        cp = f"{r['cap_per_n']:.2e}" if r['cap_per_n'] else "-"
        print(f"  {r['E_lo_eV']:.0e}-{r['E_hi_eV']:.0e} eV"
              f"{r['beam_per_pulse']:>14.3e}"
              f"{r['n_np_sim']/max(r['n_beam_sim'],1):>10.3f}"
              f"{r['rad_per_n_direct']:>14.2e}"
              f"{r['rad_per_n_npxbr']:>14.2e}"
              f"{tt:>10s}{cp:>10s}")

    if n_rad:
        print("\n--- Individual radiative-capture events ---")
        for (E, y) in sorted(acc["rad_events"]):
            print(f"  E_n = {E:11.4g} eV   cap_y = {y:8.2f} mm")


def save_results(acc, args):
    rad = np.array(acc["rad_events"], dtype=float).reshape(-1, 2)
    np.savez(args.output + ".npz",
             n_events=acc["n_events"],
             loge_edges=LOGE_EDGES,
             h_all=acc["h_all"], h_gas_np=acc["h_gas_np"],
             h_gas_rad=acc["h_gas_rad"], h_wall_al=acc["h_wall_al"],
             h_scint_h=acc["h_scint_h"],
             rad_events=rad,
             capy_edges=CAPY_EDGES,
             h_depth=np.array(acc["h_depth"]),
             budget_keys=np.array([f"{v}|{p}" for (v, p) in acc["budget"]]),
             budget_counts=np.array(list(acc["budget"].values())),
             n_per_pulse=args.n_per_pulse,
             alpha_ipc=args.alpha_ipc, br_x17=args.br_x17)
    w = args.n_per_pulse / acc["n_events"]
    summary = {
        "n_events": acc["n_events"], "n_per_pulse": args.n_per_pulse,
        "weight_per_neutron": w,
        "n_gas_np": int(acc["h_gas_np"].sum()),
        "n_gas_rad": len(acc["rad_events"]),
        "rad_events": acc["rad_events"],
        "ratio_ng_np": (len(acc["rad_events"]) / acc["h_gas_np"].sum()
                        if acc["h_gas_np"].sum() else None),
        "sigma_ratio_expected": SIGMA_RATIO,
        "rad_per_pulse_npxbr": float(acc["h_gas_np"].sum() * SIGMA_RATIO * w),
        "ipc_per_pulse": float(acc["h_gas_np"].sum() * SIGMA_RATIO * w
                               * args.alpha_ipc),
        "x17_per_pulse": float(acc["h_gas_np"].sum() * SIGMA_RATIO * w
                               * args.br_x17),
        "decades": decade_table(acc, w),
        "budget": {f"{v}|{p}": int(c) for (v, p), c in
                   acc["budget"].most_common()},
    }
    with open(args.output + ".json", "w") as jf:
        json.dump(summary, jf, indent=1)
    print(f"\nSaved: {args.output}.npz  {args.output}.json")


# ── Plot ──────────────────────────────────────────────────────────────────────
def make_plots(npz_path, out_pdf, alpha_ipc, br_x17):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    d = np.load(npz_path, allow_pickle=True)
    edges   = d["loge_edges"]
    centers = 10.0 ** ((edges[:-1] + edges[1:]) / 2)
    n_sim   = float(d["n_events"])
    w       = float(d["n_per_pulse"]) / n_sim
    h_all, h_np, h_rad = d["h_all"].astype(float), d["h_gas_np"].astype(float), \
                         d["h_gas_rad"].astype(float)

    with PdfPages(out_pdf) as pdf:
        # Page 1: beam + fates
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9), sharex=True)
        ax1.step(centers, h_all * w, where="mid", color="k")
        ax1.set_ylabel("neutrons / pulse / bin")
        ax1.set_yscale("log")
        ax1.set_title(f"Sampled EAR2 beam (sub-keV window), "
                      f"{n_sim:.2e} simulated")
        with np.errstate(divide="ignore", invalid="ignore"):
            frac_np = np.where(h_all > 0, h_np / h_all, np.nan)
        ax2.step(centers, frac_np, where="mid", label="$^3$He(n,p)t fraction")
        ax2.set_xscale("log")
        ax2.set_xlabel("$E_n$ [eV]")
        ax2.set_ylabel("fraction of beam")
        ax2.set_ylim(0, 1.05)
        ax2.legend()
        ax2.grid(alpha=0.3)
        pdf.savefig(fig); plt.close(fig)

        # Page 2: THE plot - radiative captures / pulse vs E_n
        fig, ax = plt.subplots(figsize=(9, 6.5))
        # rebin to decades for the sparse direct counts
        dec_edges = 10.0 ** np.arange(LOGE_MIN, LOGE_MAX + 1)
        dec_cent  = np.sqrt(dec_edges[:-1] * dec_edges[1:])
        rad_dec = np.array([h_rad[i*10:(i+1)*10].sum() for i in range(6)])
        np_dec  = np.array([h_np[i*10:(i+1)*10].sum() for i in range(6)])
        beam_dec = np.array([h_all[i*10:(i+1)*10].sum() for i in range(6)])

        # (n,p) x branching-ratio prediction (smooth, high stats)
        ax.step(dec_edges, np.append(np_dec, np_dec[-1]) * SIGMA_RATIO * w,
                where="post", color="C0", lw=2,
                label=r"$^3$He(n,p)t $\times\ \sigma_{n\gamma}/\sigma_{np}$"
                      " (Geant4, high stats)")
        # direct radiative counts
        has = rad_dec > 0
        ax.errorbar(dec_cent[has], rad_dec[has] * w,
                    yerr=np.sqrt(rad_dec[has]) * w, fmt="o", color="C3", ms=8,
                    capsize=4, label=r"$^3$He(n,$\gamma$) direct (Geant4)")
        for i in np.nonzero(~has)[0]:
            if beam_dec[i] > 0:
                ax.annotate("", xy=(dec_cent[i], 1.0 * w),
                            xytext=(dec_cent[i], 3.0 * w),
                            arrowprops=dict(arrowstyle="->", color="C3"))
        # Alberto thin-target overlay (per-pulse: rad/n x beam/pulse in decade)
        tt_x, tt_y, cap_x, cap_y_ = [], [], [], []
        for i, elo in enumerate(10.0 ** np.arange(LOGE_MIN, LOGE_MAX)):
            key = min(ALBERTO_THIN_TARGET, key=lambda k: abs(np.log10(k/elo))) \
                  if elo >= 1e-2 else None
            if key is not None and abs(np.log10(key / elo)) < 0.1:
                tt_x.append(dec_cent[i])
                tt_y.append(ALBERTO_THIN_TARGET[key] * beam_dec[i] * w)
                cap_x.append(dec_cent[i])
                cap_y_.append(ANALYTIC_CAP[key] * beam_dec[i] * w)
        ax.plot(tt_x, tt_y, "s--", color="C1", ms=7,
                label="thin-target table (Alberto)")
        ax.plot(cap_x, cap_y_, "^:", color="C2", ms=7,
                label=r"analytic cap $\sigma_{n\gamma}/\sigma_{np}$")

        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("$E_n$ [eV]")
        ax.set_ylabel("radiative captures / pulse / decade")
        tot = np_dec.sum() * SIGMA_RATIO * w
        ax.set_title("$^3$He(n,$\\gamma$)$^4$He per pulse — X17 production channel\n"
                     f"total: {tot:.2e} rad.capt/pulse  →  "
                     f"{tot*alpha_ipc:.2e} IPC/pulse, {tot*br_x17:.2e} X17/pulse")
        ax.legend(loc="best", fontsize=9)
        ax.grid(alpha=0.3, which="both")
        pdf.savefig(fig); plt.close(fig)

        # Page 3: absorption depth
        fig, ax = plt.subplots(figsize=(8, 5.5))
        capy_edges = d["capy_edges"]
        capy_c = (capy_edges[:-1] + capy_edges[1:]) / 2
        labels = [f"{lo:g}–{hi:g} eV" for (lo, hi) in DEPTH_SLICES]
        for i, lab in enumerate(labels):
            h = d["h_depth"][i].astype(float)
            if h.sum() > 0:
                ax.step(capy_c, h / h.sum(), where="mid", label=lab)
        ax.set_xlabel("capture y [mm]  (beam enters from −y)")
        ax.set_ylabel("normalised (n,p)t capture depth")
        ax.set_title("Self-shielding: absorption depth vs $E_n$")
        ax.legend(); ax.grid(alpha=0.3)
        pdf.savefig(fig); plt.close(fig)

    print(f"Saved: {out_pdf}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    if args.plot_only:
        make_plots(args.plot_only,
                   str(Path(args.plot_only).with_suffix("")) + ".pdf",
                   args.alpha_ipc, args.br_x17)
        return

    files = resolve_files(args.inputs, args.max_files)
    if not files:
        print("ERROR: no input files found"); sys.exit(1)
    print(f"Scanning {len(files)} files with {args.workers} workers ...")

    work = [(f, args.step_size) for f in files]
    if args.workers > 1:
        with Pool(args.workers) as pool:
            accs = pool.map(scan_file, work)
    else:
        accs = [scan_file(wk) for wk in work]
    acc = merge_accs(accs)

    print_summary(acc, args)
    save_results(acc, args)
    if args.plot:
        make_plots(args.output + ".npz", args.output + ".pdf",
                   args.alpha_ipc, args.br_x17)


if __name__ == "__main__":
    main()
