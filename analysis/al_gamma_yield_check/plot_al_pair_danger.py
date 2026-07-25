#!/usr/bin/env python3
"""Reworked two-panel 'why the pair branch is the dangerous one' schematic.

Restores the message of the 2026-07-21 leg_mechanism_schematic (84% Compton /
16% pair) but with the CORRECTED topology from the 07-23 gamma-source truth run:
the majority leg is ONE hard-Compton electron (not the old 'double Compton'),
and the minority ~18% is a genuine e+e- pair -- the same final state as the
X17/IPC signal, and mostly converting in the aluminium capsule itself.

Fractions are read from the same truth files as leg_mechanism_corrected.pdf
(docs/al_gamma_yield_check/gsrc_mechanism_nose.json, leg_mechanism_nose.json;
NOSE-FIRST 2026-07-24) so this figure stays in sync with slide 15.

Output: docs/al_gamma_yield_check/figs/al_pair_danger.pdf
"""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch

OUT = Path(__file__).resolve().parent.parent.parent / "docs/al_gamma_yield_check"
FIG = OUT / "figs"
BLUE, RED, GREEN, GRAY = "#2a78d6", "#e34948", "#008300", "#6b6a66"
plt.rcParams.update({"font.size": 10, "figure.dpi": 150})

g = json.load(open(OUT / "gsrc_mechanism_nose.json"))
n = json.load(open(OUT / "leg_mechanism_nose.json"))
vp = g["same_track_by_volproc"]
n_st = g["topology"]["same_track"]
f_same = n["same_track"] / n["legs"]


def proc_frac(proc):
    return sum(c for k, c in vp.items() if k.endswith("|" + proc)) / n_st


def loc_of(proc):
    tot = sum(c for k, c in vp.items() if k.endswith("|" + proc))
    out = {}
    for k, c in vp.items():
        v, p = k.split("|")
        if p == proc:
            out[v] = out.get(v, 0) + c / tot
    return out


fA = f_same * proc_frac("compt")    # hard Compton, single electron
fB = f_same * proc_frac("conv")     # pair production, real e+e-
locB = loc_of("conv")
al_conv = locB.get("Al capsule", 0)         # in the aluminium metal itself
alish_conv = al_conv + locB.get("capsule CFRP", 0)   # + the CFRP wrap on it


def draw_stack(ax):
    ax.add_patch(Circle((0.07, 0.60), 0.05, transform=ax.transAxes,
                        color="0.55", zorder=3))
    ax.text(0.07, 0.44, "He-3 +\nAl capsule", transform=ax.transAxes,
            ha="center", fontsize=8.5)
    for x0, wdt, c, lab in [(0.40, 0.03, "#8a6d3b", "MM+PCB"),
                            (0.60, 0.02, "#7fb3d5", "SiPM bar\n3 mm"),
                            (0.80, 0.06, "#f0b8a8", "plastic\n2 cm")]:
        ax.add_patch(Rectangle((x0, 0.30), wdt, 0.62, transform=ax.transAxes,
                               color=c, zorder=2))
        ax.text(x0 + wdt / 2, 0.24, lab, transform=ax.transAxes,
                ha="center", va="top", fontsize=8.5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")


def arrow(ax, x0, y0, x1, y1, color, ls="-", lw=2.0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), transform=ax.transAxes,
                                 arrowstyle="-|>", mutation_scale=13,
                                 color=color, lw=lw, linestyle=ls, zorder=5))


fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.6, 4.7))

# ── LEFT: the benign majority — one hard Compton electron ────────────────────
draw_stack(axL)
axL.text(0.015, 0.94, f"{fA*100:.0f}%", transform=axL.transAxes, fontsize=30,
         fontweight="bold", color=GRAY, va="top")
axL.text(0.015, 0.74, "of Al legs", transform=axL.transAxes, fontsize=10,
         color=GRAY, va="top")
axL.set_title("Compton: ONE hard e$^-$  —  a single track",
              fontsize=11, color=GRAY)
arrow(axL, 0.12, 0.60, 0.405, 0.635, BLUE)
axL.text(0.25, 0.73, "cascade $\\gamma$\n2.6–7.7 MeV", transform=axL.transAxes,
         ha="center", fontsize=8.5, color=BLUE)
axL.plot([0.415], [0.637], marker="*", ms=15, color="#eda100",
         transform=axL.transAxes, zorder=6)
arrow(axL, 0.42, 0.64, 0.95, 0.78, GREEN)
axL.text(0.72, 0.84, "e$^-$ 3–6 MeV", transform=axL.transAxes, fontsize=9,
         color=GREEN)
axL.text(0.5, 0.08, "one electron fires BOTH detectors — NOT a pair,\n"
         "so it is separable from the e$^+$e$^-$ signal",
         transform=axL.transAxes, ha="center", va="center", fontsize=8.6,
         color=GRAY, style="italic")

# ── RIGHT: the dangerous minority — a real e+e- pair, born in the Al ─────────
draw_stack(axR)
axR.text(0.015, 0.94, f"{fB*100:.0f}%", transform=axR.transAxes, fontsize=30,
         fontweight="bold", color=RED, va="top")
axR.text(0.015, 0.74, "of Al legs\n(was 16%)", transform=axR.transAxes,
         fontsize=10, color=RED, va="top")
axR.set_title("Pair production: a REAL e$^+$e$^-$  —  the dangerous one",
              fontsize=11, color=RED)
# conversion star placed AT the capsule to show it happens in the aluminium
arrow(axR, 0.12, 0.60, 0.165, 0.60, BLUE)
axR.plot([0.175], [0.60], marker="*", ms=17, color="#eda100",
         transform=axR.transAxes, zorder=6)
axR.text(0.175, 0.40, "$\\gamma\\to$ e$^+$e$^-$\nin the Al", transform=axR.transAxes,
         ha="center", fontsize=8.5, color=RED)
arrow(axR, 0.19, 0.61, 0.95, 0.82, RED)
arrow(axR, 0.19, 0.59, 0.95, 0.40, GREEN)
axR.text(0.86, 0.87, "e$^+$", transform=axR.transAxes, fontsize=11, color=RED)
axR.text(0.86, 0.34, "e$^-$", transform=axR.transAxes, fontsize=11, color=GREEN)
axR.text(0.5, 0.08, f"a genuine e$^+$e$^-$ pair = the SAME final state as X17 / IPC;\n"
         f"{al_conv*100:.0f}% converts in the Al capsule itself "
         f"({alish_conv*100:.0f}% incl. its CFRP wrap)",
         transform=axR.transAxes, ha="center", va="center", fontsize=8.6,
         color=RED, style="italic")

fig.tight_layout()
fig.savefig(FIG / "al_pair_danger.pdf")
print(f"fA(Compton)={fA:.3f}  fB(pair)={fB:.3f}  "
      f"Al-conv={al_conv:.3f}  Al+CFRP-conv={alish_conv:.3f}")
print("wrote", FIG / "al_pair_danger.pdf")
