#!/usr/bin/env python3
"""plot_particle_composition.py — what particle deposits the background energy in
the SiPM and plastic (from count_particle_composition.py). Horizontal stacked
bar of the in-gate energy-deposition fraction by particle, per detector."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent.parent / "analysis/thermal_2cm"
d = np.load(OUT / "particle_composition_2cm.npz", allow_pickle=True)
edep = d["edep"]; labels = [str(x) for x in d["labels"]]; dets = [str(x) for x in d["dets"]]

# keep particles with >=0.1% in any detector, lump the rest
frac = edep / edep.sum(axis=1, keepdims=True)
keep = np.where(frac.max(axis=0) >= 0.001)[0]
COL = {"e-": "#1f77b4", "e+": "#d62728", "proton": "#2ca02c",
       "deuteron": "#9467bd", "alpha": "#8c564b", "gamma": "#ff7f0e",
       "triton": "#e377c2", "neutron": "#7f7f7f", "other": "#bcbd22"}

fig, ax = plt.subplots(figsize=(9.5, 2.8))
y = np.arange(len(dets))
left = np.zeros(len(dets))
for j in keep:
    vals = frac[:, j] * 100
    ax.barh(y, vals, left=left, color=COL.get(labels[j], "#999"),
            label=labels[j], edgecolor="white")
    for i in range(len(dets)):
        if vals[i] >= 3:
            ax.text(left[i] + vals[i] / 2, y[i], f"{labels[j]}\n{vals[i]:.1f}%",
                    ha="center", va="center", fontsize=8,
                    color="white" if vals[i] > 8 else "black")
    left += vals
ax.set_yticks(y); ax.set_yticklabels(dets)
ax.set_xlabel("in-gate energy-deposition fraction [%]")
ax.set_xlim(0, 100)
ax.set_title("What particle deposits the background energy (sim, _2cm, in-gate)")
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, title="depositing\nparticle")
fig.text(0.5, -0.08, "e$^-$/e$^+$ = secondary electrons/positrons of the capture "
         "$\\gamma$ (Compton + pair production of the 7.72 MeV Al line). "
         "Direct proton recoil is <0.3%.", ha="center", fontsize=8.5, color="0.35")
fig.savefig(OUT / "particle_composition.pdf", bbox_inches="tight")
fig.savefig(OUT / "particle_composition.png", dpi=150, bbox_inches="tight")
print(f"wrote {OUT}/particle_composition.pdf/.png")
