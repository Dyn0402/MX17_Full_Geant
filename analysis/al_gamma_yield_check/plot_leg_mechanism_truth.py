#!/usr/bin/env python3
"""Corrected leg-mechanism schematic + truth breakdown figures.

Inputs (docs/al_gamma_yield_check/):  [NOSE-FIRST, 2026-07-24]
  gsrc_mechanism_nose.json / .npz  — gamma-source truth run (origin_vol/proc per
                                leg, 1e7 single-gamma events, 27,226 legs)
  leg_mechanism_nose.json     — neutron-mode classification (real cascade mix;
                                topology + MM-crossing fractions; 6e7 n, 1616 legs)

Outputs figs/leg_mechanism_corrected.pdf  (4-panel schematic w/ fractions)
        figs/leg_origin_breakdown.pdf     (vol x proc bars, birth-r map,
                                           per-line leg efficiency)
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "docs/al_gamma_yield_check"
FIG = OUT / "figs"
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#6b6a66"
RED, GREEN = "#e34948", "#008300"
plt.rcParams.update({"font.size": 10, "figure.dpi": 150})

g = json.load(open(OUT / "gsrc_mechanism_nose.json"))
n = json.load(open(OUT / "leg_mechanism_nose.json"))
d = np.load(OUT / "gsrc_mechanism_nose.npz")

# ── fractions ───────────────────────────────────────────────────────────────
vp = g["same_track_by_volproc"]
n_st_g = g["topology"]["same_track"]

def fsum(pred):
    return sum(c for k, c in vp.items() if pred(*k.split("|"))) / n_st_g

f_compt_g = fsum(lambda v, p: p == "compt")
f_conv_g = fsum(lambda v, p: p == "conv")

# topology from the real cascade (neutron mode)
nl = n["legs"]
f_same_n = n["same_track"] / nl
f_sgam_n = n["same_parent"] / nl
f_indep_n = n["independent"] / nl
fA = f_same_n * f_compt_g          # Compton, one crossing electron
fB = f_same_n * f_conv_g           # pair production
fC, fD = f_sgam_n, f_indep_n

def locsplit(proc):
    tot = sum(c for k, c in vp.items() if k.endswith("|" + proc))
    out = {}
    for k, c in vp.items():
        v, p = k.split("|")
        if p == proc:
            out[v] = out.get(v, 0) + c / tot
    return out

locA, locB = locsplit("compt"), locsplit("conv")
def loctxt(loc):
    order = ["Al capsule", "capsule CFRP", "MM+PCB", "SiPM bar", "air"]
    return "born in:  " + ",  ".join(
        f"{v.replace('capsule CFRP','CFRP')} {100*loc.get(v,0):.0f}%" for v in order
        if loc.get(v, 0) > 0.015)


# ── schematic ───────────────────────────────────────────────────────────────
def draw_stack(ax):
    ax.add_patch(Circle((0.06, 0.52), 0.045, transform=ax.transAxes,
                        color="0.55", zorder=3))
    ax.text(0.06, 0.38, "He-3 +\nAl capsule", transform=ax.transAxes,
            ha="center", fontsize=8)
    for x0, wdt, c, lab in [(0.36, 0.025, "#8a6d3b", "MM+PCB"),
                            (0.56, 0.018, "#7fb3d5", "SiPM bar\n3 mm"),
                            (0.76, 0.05, "#f0b8a8", "plastic\n2 cm")]:
        ax.add_patch(Rectangle((x0, 0.14), wdt, 0.74, transform=ax.transAxes,
                               color=c, zorder=2))
        ax.text(x0 + wdt / 2, 0.045, lab, transform=ax.transAxes,
                ha="center", fontsize=8)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

def arrow(ax, x0, y0, x1, y1, color, ls="-", lw=1.8, z=5):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), transform=ax.transAxes,
                                 arrowstyle="-|>", mutation_scale=11,
                                 color=color, lw=lw, linestyle=ls, zorder=z))

fig, axs = plt.subplots(2, 2, figsize=(11.4, 6.6))
panels = [
    ("A", "hard Compton, ONE e$^-$ crosses bar + plastic", fA, axs[0, 0]),
    ("B", "pair production, e$^+$/e$^-$ cross bar + plastic", fB, axs[0, 1]),
    ("C", "same $\\gamma$ interacts twice (old ``double Compton'')", fC, axs[1, 0]),
    ("D", "two different cascade $\\gamma$s, same arm", fD, axs[1, 1]),
]
for tag, title, frac, ax in panels:
    draw_stack(ax)
    ax.set_title(f"({tag})  {title}", fontsize=10.5)
    ax.text(0.985, 0.87, f"{frac*100:.1f}%", transform=ax.transAxes,
            ha="right", fontsize=21, fontweight="bold",
            color=BLUE if frac > 0.1 else GRAY)

# A: three example vertices (capsule skin / PCB / bar), one electron through
ax = axs[0, 0]
arrow(ax, 0.105, 0.53, 0.362, 0.575, BLUE)
ax.text(0.18, 0.64, "cascade $\\gamma$\n2.6--7.7 MeV", transform=ax.transAxes,
        fontsize=8, color=BLUE, ha="center")
ax.plot([0.372], [0.578], marker="*", ms=13, color="#eda100",
        transform=ax.transAxes, zorder=6)
arrow(ax, 0.375, 0.579, 0.93, 0.68, GREEN)
ax.text(0.70, 0.74, "e$^-$ 3--6 MeV", transform=ax.transAxes, fontsize=8.5,
        color=GREEN)
ax.text(0.5, 0.20, loctxt(locA), transform=ax.transAxes, fontsize=7.6,
        ha="center", color=GRAY)
# B
ax = axs[0, 1]
arrow(ax, 0.105, 0.52, 0.555, 0.52, BLUE)
ax.plot([0.567], [0.52], marker="*", ms=13, color="#eda100",
        transform=ax.transAxes, zorder=6)
arrow(ax, 0.569, 0.52, 0.93, 0.66, RED)
arrow(ax, 0.569, 0.52, 0.93, 0.40, GREEN)
ax.text(0.80, 0.71, "e$^+$", transform=ax.transAxes, fontsize=9, color=RED)
ax.text(0.80, 0.31, "e$^-$", transform=ax.transAxes, fontsize=9, color=GREEN)
ax.text(0.5, 0.20, loctxt(locB), transform=ax.transAxes, fontsize=7.6,
        ha="center", color=GRAY)
# C
ax = axs[1, 0]
arrow(ax, 0.105, 0.52, 0.555, 0.57, BLUE)
ax.plot([0.567], [0.575], marker="*", ms=12, color="#eda100",
        transform=ax.transAxes, zorder=6)
arrow(ax, 0.567, 0.575, 0.65, 0.86, GREEN, lw=1.2)
arrow(ax, 0.567, 0.575, 0.775, 0.43, BLUE, ls="--", lw=1.2)
ax.plot([0.783], [0.425], marker="*", ms=12, color="#eda100",
        transform=ax.transAxes, zorder=6)
arrow(ax, 0.785, 0.425, 0.86, 0.70, GREEN, lw=1.2)
ax.text(0.66, 0.30, "scattered $\\gamma$", transform=ax.transAxes,
        fontsize=7.6, color=BLUE)
# D
ax = axs[1, 1]
arrow(ax, 0.105, 0.55, 0.555, 0.68, BLUE)
arrow(ax, 0.105, 0.49, 0.765, 0.33, BLUE, ls="--")
ax.plot([0.567], [0.683], marker="*", ms=11, color="#eda100",
        transform=ax.transAxes, zorder=6)
ax.plot([0.773], [0.328], marker="*", ms=11, color="#eda100",
        transform=ax.transAxes, zorder=6)
arrow(ax, 0.567, 0.683, 0.63, 0.90, GREEN, lw=1.2)
arrow(ax, 0.775, 0.328, 0.85, 0.56, GREEN, lw=1.2)

fig.suptitle("How an Al capture $\\gamma$ makes a SiPM$\\wedge$plastic leg — "
             "Geant4 truth (topology: neutron run, 1616 legs; "
             "origins: $\\gamma$-source truth run, 27,226 legs)", fontsize=11.5)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(FIG / "leg_mechanism_corrected.pdf")
plt.close(fig)

# ── breakdown: volume x process | birth radius | per-line efficiency ────────
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12.6, 3.7),
                                    gridspec_kw={"width_ratios": [1.25, 1, 1]})
items = [(k, c) for k, c in vp.items() if c / n_st_g > 0.01]
other = n_st_g - sum(c for _, c in items)
labs = [k.replace("|compt", " · Compton").replace("|conv", " · pair")
        for k, _ in items] + ["everything else"]
vals = [c / n_st_g * 100 for _, c in items] + [other / n_st_g * 100]
cols = [BLUE if "compt" in k else (RED if "conv" in k else GRAY)
        for k, _ in items] + [GRAY]
y = np.arange(len(labs))[::-1]
ax1.barh(y, vals, color=cols, height=0.62)
for yy, v in zip(y, vals):
    ax1.text(v + 0.5, yy, f"{v:.1f}%", va="center", fontsize=8)
ax1.set_yticks(y, labs, fontsize=8.5)
ax1.set_xlim(0, 36)
ax1.set_xlabel("fraction of single-track legs [%]\n(birth volume · process)",
               fontsize=9)
ax1.spines[["top", "right"]].set_visible(False)

r = d["o_r"] / 10.0
ax2.hist(r, bins=np.arange(0, 48, 0.75), color=BLUE, alpha=0.85)
for x0, lab in [(1.15, "capsule"), (20.4, "MM"), (23.5, "PCB"),
                (33.2, "SiPM bar"), (41.9, "plastic")]:
    ax2.axvline(x0, color=GRAY, lw=0.6, ls=":")
    ax2.text(x0 - 0.4, ax2.get_ylim()[1] * 0.97, lab, rotation=90, va="top",
             ha="right", fontsize=7.5, color=GRAY)
ax2.set_xlabel("birth radius of the leg-making track [cm]", fontsize=9)
ax2.set_ylabel("legs")
ax2.spines[["top", "right"]].set_visible(False)

leg_E = d["lineE_leg"]
lines, ngen = d["line_E"], d["line_Ngen"]     # per-line generation counts
eff, err = [], []
for E, ng in zip(lines, ngen):
    nlg = (np.abs(leg_E - E) < 1e-3).sum()
    eff.append(nlg / ng); err.append(np.sqrt(max(nlg, 1)) / ng)
ax3.errorbar(lines, np.array(eff) * 1e3, yerr=np.array(err) * 1e3, fmt="o",
             ms=4, color=BLUE, lw=1)
ax3.axvline(7.724, color=GRAY, lw=0.6, ls=":")
ax3.text(7.6, 0.5, "$\\gamma_0$ 7.724:\n49% of legs", fontsize=8, ha="right",
         color=GRAY)
ax3.set_xlabel("cascade line energy [MeV]", fontsize=9)
ax3.set_ylabel("legs per emitted $\\gamma$  [$\\times10^{-3}$]")
ax3.set_ylim(0, 6.2)
ax3.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIG / "leg_origin_breakdown.pdf")
print("A/B/C/D fractions: %.3f %.3f %.3f %.3f" % (fA, fB, fC, fD))
print("locA:", {k: round(v, 3) for k, v in locA.items()})
print("locB:", {k: round(v, 3) for k, v in locB.items()})
print("wrote figs")
