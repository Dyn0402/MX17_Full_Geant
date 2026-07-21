#!/usr/bin/env python3
"""plot_event_display.py — 2D event display of a real Al-capture leg event from a
trajdump CSV. Top-down (x-z) and side (x-y) projections; trajectory segments
colored by particle, geometry bands overlaid. Shows how the 7.72 MeV Al capture
gamma deposits in both the SiPM wall and the plastic of one arm."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parent.parent / "analysis/thermal_2cm"
rows = list(csv.reader(open(OUT / "al_leg_event.csv")))[1:]

COL = {"gamma": "#1f77b4", "e-": "#2ca02c", "e+": "#d62728",
       "neutron": "#999999", "proton": "#8c564b"}
LW  = {"gamma": 1.3, "neutron": 1.3}

def draw(ax, ia, ib, xlabel, ylabel):
    H = 260
    ax.add_patch(Circle((0, 0), 13, color="0.6", zorder=1))                              # capsule
    ax.add_patch(Rectangle((-234, -H), 30, 2*H, color="#8fbf8f", alpha=0.30, zorder=0))  # MM drift gas (active)
    ax.add_patch(Rectangle((-240, -H), 6, 2*H, color="#8a6d3b", alpha=0.65, zorder=0))   # MM PCB (solid scatterer)
    ax.add_patch(Rectangle((-349, -H), 35, 2*H, color="#7fb3d5", alpha=0.12, zorder=0))  # SiPM container
    ax.add_patch(Rectangle((-333, -H), 3, 2*H, color="#1f77b4", alpha=0.75, zorder=0))   # SiPM active 3 mm
    ax.add_patch(Rectangle((-430, -H), 20, 2*H, color="#e6866f", alpha=0.45, zorder=0))  # plastic 2 cm
    for r in rows:
        p = r[3]
        ax.plot([float(r[ia]), float(r[ia + 3])], [float(r[ib]), float(r[ib + 3])],
                color=COL.get(p, "k"), lw=LW.get(p, 0.4),
                alpha=0.65 if p in ("e-", "e+") else 0.9, zorder=2 if p in ("e-", "e+") else 3)
    # conv (pair-production) vertex
    for r in rows:
        if r[14] == "conv":
            ax.plot(float(r[ia]), float(r[ib]), "*", ms=16, color="gold",
                    mec="k", mew=1, zorder=5)
            break
    yl = H + 16
    ax.text(-331, yl, "SiPM\n3 mm", ha="center", fontsize=7.5, color="#1f5f8b")
    ax.text(-420, yl, "plastic\n2 cm", ha="center", fontsize=7.5, color="#a03a26")
    ax.text(-219, yl, "MM drift\n30 mm", ha="center", fontsize=7, color="#3a6b3a")
    ax.text(-237, -yl, "MM PCB\n(Cu/FR4)", ha="center", fontsize=7, color="#8a6d3b", va="top")
    ax.text(0, 32, "He-3 +\nAl capsule", ha="center", fontsize=7.5)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_xlim(-460, 70); ax.set_ylim(-H-58, H+58); ax.grid(True, alpha=0.2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.2))
draw(ax1, 5, 7, "x  [mm]   (−X arm ←)", "z  [mm]")     # top-down x-z
ax1.set_title("Top-down (x–z)")
draw(ax2, 5, 6, "x  [mm]   (−X arm ←)", "y  [mm]  (beam)")  # side x-y
ax2.set_title("Side (x–y, beam vertical)")

legend = [Line2D([0], [0], color=COL["gamma"], lw=2, label="γ"),
          Line2D([0], [0], color=COL["e-"], lw=2, label="e−"),
          Line2D([0], [0], color=COL["e+"], lw=2, label="e+"),
          Line2D([0], [0], color=COL["neutron"], lw=2, label="neutron"),
          Line2D([0], [0], marker="*", color="gold", mec="k", lw=0, ms=12,
                 label="pair-production (conv)")]
fig.legend(handles=legend, loc="upper center", ncol=5, fontsize=9,
           bbox_to_anchor=(0.5, 1.02))
fig.suptitle("Real Geant4 event: $^{27}$Al(n,γ) 7.72 MeV γ → pair production → "
             "EM shower deposits in SiPM (1.88 MeV) + plastic (2.70 MeV) of the −X arm",
             y=1.07, fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "event_display_al_leg.pdf", bbox_inches="tight")
fig.savefig(OUT / "event_display_al_leg.png", dpi=145, bbox_inches="tight")
print(f"wrote {OUT}/event_display_al_leg.pdf/.png")
