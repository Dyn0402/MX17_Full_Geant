#!/usr/bin/env python3
"""
plot_timedist_bysource_thermal.py — stacked-by-source time distributions.
Reads count_timedist_by_source_thermal.py npz and makes a 3-panel figure
(SiPM-wall singles / plastic singles / legs) with a shared LINEAR time axis;
within each panel the rate is a stacked area broken down by neutron-capture
source. Also writes a compact per-source CSV for overlay.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
NPZ = Path(sys.argv[1]) if len(sys.argv) > 1 else \
      ROOT / "analysis/thermal_2cm/timedist_bysource_2cm.npz"
OUT = ROOT / "analysis/thermal_2cm"

d = np.load(NPZ, allow_pickle=True)
tedges = d["tedges"]; tc = d["tc"]; dt = np.diff(tedges)
n_ev = float(d["n_events"]); n_pp = float(d["n_pulse"]); w = n_pp / n_ev
sources = [str(s) for s in d["sources"]]
COLORS = ["#4C72B0", "#8C8C8C", "#DD8452", "#C44E52", "#55A868", "#DADADA"]

panels = [("SiPM-wall singles", "sipm"),
          ("Plastic singles",   "plas"),
          ("Legs (SiPM $\\wedge$ plastic / arm)", "leg")]

# rate per source: (nsrc, ntbin) -> per pulse per ms
def rates(key):
    return d[key].astype(float) * w / dt        # (nsrc, ntbin)

# order sources by total contribution (largest at bottom of stack), by legs
order = np.argsort(-d["leg"].astype(float).sum(axis=1))
tot_sipm = rates("sipm").sum(0)
tmax = float(tc[tot_sipm > 1e-3 * tot_sipm.max()].max()) if tot_sipm.max() > 0 else tc.max()
mt = tc <= tmax

fig, axes = plt.subplots(3, 1, figsize=(9.0, 9.4), sharex=True)
for ax, (label, key) in zip(axes, panels):
    R = rates(key)                                    # (nsrc, ntbin)
    stack = [R[s][mt] for s in order]
    labels = [sources[s] for s in order]
    cols   = [COLORS[s % len(COLORS)] for s in order]
    ax.stackplot(tc[mt], *stack, labels=labels, colors=cols,
                 edgecolor="none", alpha=0.9)
    ax.set_ylabel(f"{label}\n[/pulse/ms]", fontsize=9)
    ax.grid(True, alpha=0.25)
    ax.margins(y=0)
    ax.set_ylim(bottom=0)
    tot = float(d[key].sum()) * w                    # integral = counts * w
    al  = float(d[key][0].sum()) * w                 # Al-capsule contribution
    ax.text(0.015, 0.93, f"$\\int$ = {tot:,.0f}/pulse   (Al {al/tot*100:.0f}%)",
            transform=ax.transAxes, ha="left", va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.85))
axes[0].set_title("Thermal-gate rates vs arrival time, by capture source "
                  "(sim, _2cm, 0.5 MIP)", fontsize=11)
axes[0].legend(loc="upper right", fontsize=8, ncol=2, framealpha=0.9)
axes[-1].set_xlabel("neutron arrival time  t [ms]   (TOF, 19.5 m; gate t > 1 ms)")
axes[-1].set_xlim(1.0, tmax)
fig.align_ylabels(axes)
fig.tight_layout()
fig.savefig(OUT / "timedist_bysource.pdf", bbox_inches="tight")
fig.savefig(OUT / "timedist_bysource.png", dpi=140, bbox_inches="tight")

# compact CSV: one block per quantity, columns = sources
with open(OUT / "timedist_bysource.csv", "w") as f:
    f.write("# thermal-gate sim rates vs arrival time by capture source, "
            "per pulse per ms, 4 walls, 0.5 MIP, plastic@2.0cm. "
            "n_events=%d n_pulse=%.4e\n" % (n_ev, n_pp))
    hdr = ["quantity", "t_lo_ms", "t_hi_ms", "t_mid_ms"] + sources + ["total"]
    f.write(",".join(hdr) + "\n")
    for label, key in panels:
        R = rates(key)
        for i in range(len(tc)):
            row = [key, f"{tedges[i]:.4e}", f"{tedges[i+1]:.4e}", f"{tc[i]:.4e}"]
            row += [f"{R[s][i]:.4e}" for s in range(len(sources))]
            row += [f"{R[:,i].sum():.4e}"]
            f.write(",".join(row) + "\n")

print("legs/pulse by source:")
for s in order:
    print(f"  {sources[s]:16s} {d['leg'][s].sum()*w:7.1f}")
print(f"wrote {OUT}/timedist_bysource.pdf/.png/.csv")
