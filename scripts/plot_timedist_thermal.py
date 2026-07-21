#!/usr/bin/env python3
"""
plot_timedist_thermal.py — stacked time distributions for data overlay
======================================================================
Reads the histogram npz from count_timedist_thermal.py and produces:
  - a stacked, shared-time-axis figure (SiPM-wall singles / plastic singles /
    per-arm SiPM^plastic legs), each on its own vertical axis, in-gate (t>1 ms)
  - a compact CSV (per-pulse-per-ms rates + Poisson errors) for overlaying the
    measured rates.

Rates are per beam pulse, summed over the 4 walls, at 0.5 MIP.  Plastic/leg
lines are drawn for the 2.0 cm (current) MIP calibration; the 2.5 cm value is
included in the CSV as a bracket on the plastic-threshold energy.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NPZ  = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "analysis/thermal_2cm/timedist_2cm.npz"
OUT  = ROOT / "analysis/thermal_2cm"

d = np.load(NPZ)
tedges = d["tedges"]; tc = d["tc"]; dt = np.diff(tedges)
n_ev = float(d["n_events"]); n_pp = float(d["n_pulse"]); w = n_pp / n_ev

def rate(counts):           # per pulse per ms
    return counts * w / dt
def err(counts):            # Poisson per pulse per ms
    return np.sqrt(counts) * w / dt

series = [
    ("SiPM-wall singles",           "sipm",   "#1f77b4"),
    ("Plastic singles",             "plas20", "#d62728"),
    ("Legs (SiPM $\\wedge$ plastic / arm)", "leg20", "#2ca02c"),
]

# trim dead low-energy tail (where total capture rate is negligible)
tot = d["sipm"].astype(float)
tmax = float(tc[tot > 1e-3 * tot.max()].max()) if tot.max() > 0 else tc.max()

fig, axes = plt.subplots(3, 1, figsize=(8.6, 9.2), sharex=True)
for ax, (label, key, col) in zip(axes, series):
    c = d[key].astype(float); r = rate(c); e = err(c)
    m = (tc <= tmax)
    ax.errorbar(tc[m], r[m], yerr=e[m], fmt="o-", ms=3, lw=1.4, color=col,
                capsize=0, elinewidth=0.8)
    ax.set_ylim(bottom=0)          # linear (per request)
    ax.set_ylabel(f"{label}\n[/pulse/ms]", fontsize=9)
    ax.grid(True, which="both", alpha=0.25)
    ax.text(0.985, 0.9, f"$\\int$ = {c.sum()*w:,.0f}/pulse",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec=col, alpha=0.8))
axes[0].set_title("Thermal-gate SiPM / plastic / coincidence rates vs arrival time "
                  "(sim, _2cm, 0.5 MIP)", fontsize=11)
axes[-1].set_xlabel("neutron arrival time  t [ms]   (TOF, 19.5 m; gate t > 1 ms)")
axes[-1].set_xlim(1.0, tmax)
fig.align_ylabels(axes)
fig.tight_layout()
fig.savefig(OUT / "timedist_thermal.pdf", bbox_inches="tight")
fig.savefig(OUT / "timedist_thermal.png", dpi=140, bbox_inches="tight")

# compact CSV for overlay
with open(OUT / "timedist_thermal.csv", "w") as f:
    f.write("# thermal-gate sim rates vs arrival time, per pulse, summed over 4 walls, "
            "0.5 MIP. n_events=%d n_pulse=%.4e. Poisson errors.\n" % (n_ev, n_pp))
    f.write("t_lo_ms,t_hi_ms,t_mid_ms,"
            "sipm_per_pulse_per_ms,sipm_err,"
            "plastic2.0cm_per_pulse_per_ms,plastic2.0cm_err,"
            "plastic2.5cm_per_pulse_per_ms,plastic2.5cm_err,"
            "leg2.0cm_per_pulse_per_ms,leg2.0cm_err,"
            "leg2.5cm_per_pulse_per_ms,leg2.5cm_err\n")
    for i in range(len(tc)):
        f.write(",".join(f"{x:.6e}" for x in [
            tedges[i], tedges[i+1], tc[i],
            rate(d["sipm"])[i],   err(d["sipm"])[i],
            rate(d["plas20"])[i], err(d["plas20"])[i],
            rate(d["plas25"])[i], err(d["plas25"])[i],
            rate(d["leg20"])[i],  err(d["leg20"])[i],
            rate(d["leg25"])[i],  err(d["leg25"])[i],
        ]) + "\n")

print("totals /pulse:  SiPM %.0f  plastic(2.0) %.0f  legs(2.0) %.0f"
      % (d["sipm"].sum()*w, d["plas20"].sum()*w, d["leg20"].sum()*w))
print(f"wrote {OUT}/timedist_thermal.pdf/.png/.csv")
