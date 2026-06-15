#!/usr/bin/env python3
"""
make_he4_levels_fig.py

Level diagram of 4He near the n+3He / p+3H thresholds, drawn to make the
formation-vs-decay story concrete:
  - which states exist (energies, widths, Jpi, isospin T) -- TUNL A=4
    (Tilley/Weller/Hale, NPA 541 (1992) 1);
  - the particle thresholds the captured neutron sits at;
  - which states are reached by s-wave vs p-wave from n+3He;
  - that there is NO low-lying 1+ level: the lowest 1+ is at ~28.3 MeV, so the
    M1 radiative capture goes through the non-resonant 3S1 (1+) CONTINUUM at
    threshold (direct capture), not through a 1+ resonance;
  - isospin: the ground state is T=0, so E1 capture (isovector) to it is
    isospin-forbidden except from the T=1 1- states (the GDR) -> E1 turns on
    only toward MeV.

The y-axis is BROKEN between 4 and 16 MeV: 4He has exactly one bound state (the
0+ g.s. at 0) and a 20 MeV gap up to the first excited state, so the break lets
the excited-state cluster (20-28 MeV) spread out while still showing the g.s.

Widths are shown as shaded bands -- these states are BROAD and overlapping
(4He has no sharp levels above threshold; they are continuum resonances).

Output: docs/e0_branch/figs/fig_he4_levels.{pdf,png}
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parent.parent / "docs" / "e0_branch" / "figs"

# (Ex MeV, Jpi, T, width MeV, entrance wave from n+3He, g.s. transition note)
# TUNL A=4 (Tilley, Weller, Hale, NPA 541 (1992) 1). The states above ~20 MeV
# are broad continuum resonances from the evaluated 4He R-matrix.
LEVELS = [
    (0.00,  "0^+", 0, 0.00, "",  "ground state"),
    (20.21, "0^+", 0, 0.50, "s", "E0"),
    (21.01, "0^-", 0, 0.84, "p", "forbidden"),
    (21.84, "2^-", 0, 2.01, "p", "M2 (slow)"),
    (23.33, "2^-", 1, 5.01, "p", "M2 (slow)"),
    (24.25, "1^-", 0, 6.10, "p", "E1 (T=0: hindered)"),
    (25.95, "1^-", 1, 12.66, "p", "E1 (T=1: GDR)"),
    (28.31, "1^+", 0, 9.00, "s", "M1"),
]
THRESHOLDS = [
    (19.815, r"$p+{}^{3}$H"),
    (20.578, r"$n+{}^{3}$He"),
    (23.847, r"$d+d$"),
]
# colour by role for the e+e- pair signal (transition to the 0+ ground state)
ROLE_COLOR = {"E0": "#c0392b", "M1": "#2471a3",
              "E1 (T=0: hindered)": "#16a085", "E1 (T=1: GDR)": "#16a085",
              "forbidden": "0.6", "M2 (slow)": "0.6", "ground state": "k"}
ROLE_TAG = {"E0": "s$\\to$E0", "M1": "s$\\to$M1",
            "E1 (T=0: hindered)": "p$\\to$E1 (T=0, hind.)",
            "E1 (T=1: GDR)": "p$\\to$E1 (T=1, GDR)",
            "forbidden": "p, g.s. forbidden", "M2 (slow)": "p$\\to$M2 (slow)"}

# broken axis: bottom shows the bound g.s., top shows the continuum cluster
LO_LIM = (-1.6, 4.0)
HI_LIM = (16.0, 30.5)
xL, xR = 0.18, 0.74


def draw_level(ax, Ex, jp, T, w, role):
    c = ROLE_COLOR[role]
    if w > 0:
        ax.add_patch(Rectangle((xL, Ex - w / 2), xR - xL, w,
                     facecolor=c, alpha=0.11, edgecolor="none", zorder=1))
    ax.plot([xL, xR], [Ex, Ex], color=c, lw=2.6, zorder=3)
    if Ex > 0:
        ax.text(xL - 0.015, Ex, f"{Ex:.2f}", va="center", ha="right",
                fontsize=9, color="0.3")
        ax.text(xR + 0.015, Ex, f"${jp},\\,T={T}$", va="center", ha="left",
                fontsize=10.5, color=c, fontweight="bold")
        ax.text(xR + 0.205, Ex, f"[{ROLE_TAG[role]}]", va="center",
                ha="left", fontsize=7.8, color=c, style="italic")
    else:
        ax.text(xR + 0.015, Ex, f"${jp},\\,T={T}$  (ground state)", va="center",
                ha="left", fontsize=10.5, color=c, fontweight="bold")


def main():
    fig, (axH, axL) = plt.subplots(2, 1, figsize=(9.8, 9.4), sharex=True,
                                   gridspec_kw={"height_ratios": [14.5, 5.6],
                                                "hspace": 0.06})
    axH.set_ylim(*HI_LIM)
    axL.set_ylim(*LO_LIM)

    for Ex, jp, T, w, wave, role in LEVELS:
        ax = axH if Ex >= 16 else axL
        draw_level(ax, Ex, jp, T, w, role)

    # thresholds (in the top panel)
    for Eth, lab in THRESHOLDS:
        axH.axhline(Eth, xmin=0.16, xmax=0.74, color="darkorange", ls="--",
                    lw=1.2, zorder=2)
        axH.text(0.02, Eth, f"{lab} ({Eth:.2f})", va="center",
                 ha="left", fontsize=8.2, color="darkorange")

    # teaching annotations (top panel)
    axH.annotate("lowest $1^+$ is here ($\\sim$28 MeV): $\\sim$7.7 MeV ABOVE\n"
                 "threshold, very broad. $\\Rightarrow$ the M1 capture is DIRECT,\n"
                 "through the non-resonant $^3S_1$ continuum at\n"
                 "threshold, NOT through a $1^+$ level.",
                 xy=(xR + 0.05, 28.31), xytext=(0.43, 17.4),
                 fontsize=8.2, color="#2471a3", ha="left",
                 arrowprops=dict(arrowstyle="->", color="#2471a3"))
    axH.annotate("the ONLY sub-threshold state ($-0.37$ MeV).\n"
                 "Its tail boosts E0 capture at low $E_n$.",
                 xy=(xL + 0.05, 20.21), xytext=(0.30, 16.7),
                 fontsize=8.2, color="#c0392b", ha="left",
                 arrowprops=dict(arrowstyle="->", color="#c0392b"))
    axL.annotate("the only BOUND state of $^4$He;\n"
                 "everything above is a continuum resonance",
                 xy=(xR, 0.0), xytext=(0.30, 2.3),
                 fontsize=8.2, color="0.25", ha="left",
                 arrowprops=dict(arrowstyle="->", color="0.4"))

    # --- broken-axis cosmetics ------------------------------------------------
    axH.spines["bottom"].set_visible(False)
    axL.spines["top"].set_visible(False)
    axH.tick_params(bottom=False)
    d = .008
    kw = dict(transform=axH.transAxes, color="k", clip_on=False, lw=1)
    axH.plot((-d, +d), (-d, +d), **kw); axH.plot((1 - d, 1 + d), (-d, +d), **kw)
    kw.update(transform=axL.transAxes)
    sc = 14.5 / 5.6   # scale break-mark slope to the two panel heights
    axL.plot((-d, +d), (1 - d * sc, 1 + d * sc), **kw)
    axL.plot((1 - d, 1 + d), (1 - d * sc, 1 + d * sc), **kw)

    for ax in (axH, axL):
        ax.set_xlim(0.0, 1.30)
        ax.set_xticks([])
    axL.set_xlabel("(y-axis broken 4$\\to$16 MeV: $^4$He has a 20 MeV gap "
                   "between its bound g.s. and the first excited state)",
                   fontsize=8.5, color="0.4")
    fig.text(0.045, 0.55, "excitation energy  $E_x$  [MeV]", rotation=90,
             va="center", fontsize=11)
    axH.set_title("$^4$He levels near the $n+{}^3$He threshold (TUNL A=4) --- "
                  "bands = widths, colour = how the state\ndecays to the "
                  "$0^+$ ground state. States above threshold are broad "
                  "continuum resonances.",
                  fontsize=11)

    leg = [Line2D([], [], color="#c0392b", lw=3, label="E0 $\\to$ 100% pairs ($\\gamma$-dark)"),
           Line2D([], [], color="#2471a3", lw=3, label="M1 $\\to$ $\\gamma$ + 0.35% pairs"),
           Line2D([], [], color="#16a085", lw=3, label="E1 $\\to$ $\\gamma$ + 0.35% pairs (T=1 only)"),
           Line2D([], [], color="0.6", lw=3, label="forbidden / slow (spectator)"),
           Line2D([], [], color="darkorange", ls="--", lw=1.2, label="particle threshold")]
    axL.legend(handles=leg, loc="lower right", fontsize=8.2, framealpha=0.95,
               title="g.s. transition type", title_fontsize=8.5)

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_he4_levels.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved -> {OUT}/fig_he4_levels.[pdf,png]")


if __name__ == "__main__":
    main()
