#!/usr/bin/env python3
"""Reconciliation ladder: Alberto's single-arm coincidence estimate (~9/pulse)
vs the sim's Al-attributable legs (199/pulse). Steps apply, in order:
capture/intensity truth, full 4-arm aperture, all cascade lines >2.5 MeV,
and the measured single-electron punch-through mechanism.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent.parent
FIG = ROOT / "docs/al_gamma_yield_check/figs"
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#6b6a66"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150})

rows = [
    ("Alberto: $\\gamma_0\\times$ 4.5% $\\Omega$ $\\times$ 1% Compton", 8.8, GRAY, ""),
    ("true captures + 21.3% line intensity", 3.2, ORANGE, "$\\times$0.36"),
    ("full 4-arm plastic aperture (18.5%)", 13.0, GRAY, "$\\times$4.1"),
    ("measured response to a 7.72 MeV $\\gamma$ (0.54%/$\\gamma$)", 37.8, GRAY, "$\\times$2.9"),
    ("+ other cascade lines (2.6--7.7 MeV, $\\gamma$-source run)", 105.7, GRAY, "$\\times$2.8"),
    ("+ full G4NDL cascade (multiplicity + co-adding)", 199., BLUE, "$\\times$1.9"),
    ("Geant4 legs, Al$-$noAl (truth)", 198.8, AQUA, ""),
]
fig, ax = plt.subplots(figsize=(8.8, 3.5))
ypos = np.arange(len(rows))[::-1]
for y, (lab, v, c, f) in zip(ypos, rows):
    ax.barh(y, v, color=c, height=0.62)
    ax.text(v * 1.09, y, f"{v:,.1f}" + (f"   {f}" if f else ""),
            va="center", fontsize=9)
    ax.text(0.55, y, lab, va="center", ha="right", fontsize=9.5)
ax.set_xscale("log"); ax.set_xlim(0.6, 1.5e3)
ax.set_yticks([])
ax.set_xlabel("single-arm SiPM$\\,\\wedge\\,$plastic coincidences / pulse  (in-gate)")
ax.spines["left"].set_visible(False)
fig.tight_layout()
fig.savefig(FIG / "coincidence_ladder.pdf")
print("wrote", FIG / "coincidence_ladder.pdf")
