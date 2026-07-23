#!/usr/bin/env python3
"""plot_al_noal_gcut_compare.py — two cross-check figures requested by
collaborators, built on the same "thermal-gate rate vs arrival time" metric
as analysis/thermal_2cm/timedist_thermal.png (leg = SiPM-wall AND plastic
>=0.5 MIP, per arm, summed over 4 arms, per pulse per ms, t>1 ms gate):

  1. al_vs_noal.png   — baseline (Al capsule) vs --no-al (vacuum vessel).
                         Tests whether the Al(n,g) 7.72 MeV background really
                         drives the observed in-gate rate.
  2. gamma_cut_scan.png — baseline (100 um) vs gamma production-cut 50/20/10/5 um.
                         Tests whether the trigger rate depends on the gamma
                         production-cut choice (transport-fidelity sanity check).

Reads the npz files written by count_timedist_thermal_compare.py.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
IN   = ROOT / "analysis/al_pair_crosscheck"
OUT  = IN

def load(name):
    d = np.load(IN / f"timedist_{name}.npz")
    n_ev, n_pp = float(d["n_events"]), float(d["n_pulse"])
    w = n_pp / n_ev
    tc, dt = d["tc"], np.diff(d["tedges"])
    leg = d["leg20"].astype(float)
    rate = leg * w / dt
    err  = np.sqrt(leg) * w / dt
    total_per_pulse = leg.sum() * w
    return tc, rate, err, total_per_pulse, n_ev


def make_plot(runs, out_stem, title):
    """runs: list of (label, npz_name, color)."""
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    totals = []
    for label, name, col in runs:
        tc, rate, err, tot, n_ev = load(name)
        m = rate > 0
        ax.errorbar(tc[m], rate[m], yerr=err[m], fmt="o-", ms=3, lw=1.3,
                    color=col, capsize=0, elinewidth=0.7,
                    label=f"{label}  ({tot:,.0f}/pulse, N={n_ev:,.0f})")
        totals.append((label, tot))
    ax.set_yscale("log")
    ax.set_xlabel("neutron arrival time  t [ms]   (TOF, 19.5 m; gate t > 1 ms)")
    ax.set_ylabel("legs (SiPM $\\wedge$ plastic / arm, 0.5 MIP)\n[/pulse/ms]")
    ax.set_title(title, fontsize=11)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / f"{out_stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{out_stem}.png", dpi=140, bbox_inches="tight")
    print(f"wrote {out_stem}.png/.pdf")
    print(f"  totals/pulse: " + ", ".join(f"{l}={t:.1f}" for l, t in totals))
    return totals


if __name__ == "__main__":
    make_plot(
        [("Al capsule (baseline)", "baseline_100um", "#d62728"),
         ("No Al (vacuum vessel)", "noAl",           "#1f77b4")],
        "al_vs_noal",
        "Thermal-gate background: Al capsule vs. no-Al cross-check",
    )
    make_plot(
        [("100 $\\mu$m (baseline)", "baseline_100um", "#1f77b4"),
         ("50 $\\mu$m",             "gcut50um",       "#2ca02c"),
         ("20 $\\mu$m",             "gcut20um",       "#ff7f0e"),
         ("10 $\\mu$m",             "gcut10um",       "#d62728"),
         ("5 $\\mu$m",              "gcut5um",        "#9467bd")],
        "gamma_cut_scan",
        "Thermal-gate background: sensitivity to the gamma production cut",
    )
