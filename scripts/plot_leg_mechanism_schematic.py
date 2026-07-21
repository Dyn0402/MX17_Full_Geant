#!/usr/bin/env python3
"""Schematic of the two ways the 7.72 MeV Al(n,gamma) capture gamma makes a
SiPM-and-plastic leg: (1) double Compton (typical), (2) pair-production shower
(high-energy). Complements the real event display."""
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch

OUT = Path(__file__).resolve().parent.parent / "analysis/thermal_2cm"

def arm(ax):
    ax.add_patch(Circle((0, 0), 0.5, color="0.6"))
    ax.text(0, -0.9, "He-3 +\nAl capsule", ha="center", fontsize=8)
    ax.add_patch(Rectangle((3.3, -2), 0.5, 4, color="#7fb3d5", alpha=0.5))
    ax.text(3.55, 2.35, "SiPM wall", ha="center", fontsize=8, color="#1f5f8b")
    ax.add_patch(Rectangle((5.3, -2), 0.35, 4, color="#e6866f", alpha=0.5))
    ax.text(5.48, 2.35, "plastic", ha="center", fontsize=8, color="#a03a26")
    ax.set_xlim(-1.2, 7); ax.set_ylim(-3, 3.2); ax.axis("off")

def g(ax, x0, y0, x1, y1, **k):
    ax.annotate("", (x1, y1), (x0, y0),
                arrowprops=dict(arrowstyle="-|>", color="#1f77b4", lw=2, **k))
def e(ax, x0, y0, x1, y1, col):
    ax.annotate("", (x1, y1), (x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=col, lw=1.6))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 5.2))

# (1) double Compton
arm(a1)
g(a1, 0.4, 0.1, 3.35, 0.4)                       # gamma to SiPM
e(a1, 3.55, 0.4, 3.75, 1.4, "#2ca02c")           # recoil e- in SiPM
g(a1, 3.55, 0.4, 5.35, -0.5, ls="--")            # scattered gamma to plastic
e(a1, 5.45, -0.5, 5.65, 0.6, "#2ca02c")          # recoil e- in plastic
a1.text(1.7, 0.75, "7.72 MeV γ", color="#1f77b4", fontsize=9)
a1.text(3.9, 1.4, "Compton e−", color="#2ca02c", fontsize=8)
a1.text(4.3, -0.95, "scattered γ", color="#1f77b4", fontsize=8)
a1.text(5.75, 0.7, "Compton e−", color="#2ca02c", fontsize=8)
a1.text(-1.1, 2.7, "84%", fontsize=34, fontweight="bold", color="#1f5f8b",
        va="top", ha="left")
a1.text(0.7, 2.55, "of Al legs", fontsize=10, color="#1f5f8b", va="top")
a1.set_title("(1) Double Compton  —  the typical leg\n"
             "γ Compton-scatters in SiPM, the scattered γ Compton-scatters in plastic",
             fontsize=10)

# (2) pair production + shower
arm(a2)
g(a2, 0.4, 0.1, 3.3, 0.15)                        # gamma into SiPM region
a2.plot(3.35, 0.15, "*", ms=17, color="gold", mec="k", zorder=5)
e(a2, 3.4, 0.2, 4.6, 1.6, "#d62728")              # e+ 
e(a2, 3.4, 0.1, 4.7, -1.4, "#2ca02c")             # e-
e(a2, 4.6, 1.6, 5.5, 1.1, "#d62728")              # e+ into plastic
e(a2, 4.7, -1.4, 5.4, -0.4, "#2ca02c")            # e- into plastic
g(a2, 4.6, 1.6, 6.6, 2.2, ls=":")                 # brems gamma escaping
g(a2, 4.7, -1.4, 6.6, -2.1, ls=":")
a2.text(1.7, 0.5, "7.72 MeV γ", color="#1f77b4", fontsize=9)
a2.text(3.0, 0.9, "pair prod.\n(conv)", fontsize=8, ha="center")
a2.text(4.75, 1.75, "e+", color="#d62728", fontsize=9)
a2.text(4.85, -1.65, "e−", color="#2ca02c", fontsize=9)
a2.text(6.2, 2.35, "brems γ", color="#1f77b4", fontsize=8)
a2.text(-1.1, 2.7, "16%", fontsize=34, fontweight="bold", color="#a03a26",
        va="top", ha="left")
a2.text(0.7, 2.55, "of Al legs\n(higher energy)", fontsize=10, color="#a03a26", va="top")
a2.set_title("(2) Pair production + EM shower  —  the high-energy leg\n"
             "7.72 MeV γ converts (mostly in the Al capsule); e± spray both detectors",
             fontsize=10)

fig.suptitle("How the $^{27}$Al(n,γ) 7.72 MeV capture γ makes a SiPM∧plastic leg "
             "(both deposit >0.5 MIP in one arm)", y=1.02, fontsize=12)
fig.text(0.5, -0.02, "Detector interactions are Compton-dominated (thin low-Z); "
         "pair production happens mainly in the denser Al capsule and dominates the "
         "highest-energy legs. The soft 2.22 MeV H-capture γ rarely clears both thresholds.",
         ha="center", fontsize=8.5, color="0.35")
fig.tight_layout()
fig.savefig(OUT / "leg_mechanism_schematic.pdf", bbox_inches="tight")
fig.savefig(OUT / "leg_mechanism_schematic.png", dpi=145, bbox_inches="tight")
print(f"wrote {OUT}/leg_mechanism_schematic.pdf/.png")
