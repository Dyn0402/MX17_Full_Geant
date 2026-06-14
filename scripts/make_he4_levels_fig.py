#!/usr/bin/env python3
"""
make_he4_levels_fig.py

Level diagram of 4He near the n+3He / p+3H thresholds, drawn to make the
formation-vs-decay story concrete:
  - which states exist (energies, widths, Jpi) -- TUNL A=4 (Tilley/Weller/Hale,
    NPA 541 (1992) 1);
  - the thresholds the captured neutron sits at;
  - which states are reached by s-wave vs p-wave from n+3He;
  - and, crucially, that there is NO low-lying 1+ level: the lowest 1+ is at
    ~28.3 MeV.  The M1 radiative capture therefore goes through the
    non-resonant 3S1 (1+) CONTINUUM at threshold (direct capture), not through
    a 1+ resonance.

Widths are shown as shaded bands -- these states are BROAD and overlapping,
which is itself part of the story (4He has no sharp levels above threshold).

Output: docs/e0_branch/figs/fig_he4_levels.{pdf,png}
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = Path(__file__).resolve().parent.parent / "docs" / "e0_branch" / "figs"

# (Ex MeV, Jpi, T, width MeV, entrance wave from n+3He, g.s. transition note)
LEVELS = [
    (0.00,  "0^+", 0, 0.00, "",  "ground state"),
    (20.21, "0^+", 0, 0.50, "s", "E0"),
    (21.01, "0^-", 0, 0.84, "p", "forbidden"),
    (21.84, "2^-", 0, 2.01, "p", "M2 (slow)"),
    (23.33, "2^-", 1, 5.01, "p", "M2 (slow)"),
    (24.25, "1^-", 0, 6.10, "p", "E1"),
    (25.95, "1^-", 1, 12.66, "p", "E1"),
    (28.31, "1^+", 0, 9.00, "s", "M1"),
]
THRESHOLDS = [
    (19.815, r"$p+{}^{3}$H"),
    (20.578, r"$n+{}^{3}$He"),
    (23.847, r"$d+d$"),
]
# colour by role for the e+e- pair signal (transition to the 0+ ground state)
ROLE_COLOR = {"E0": "#c0392b", "M1": "#2471a3", "E1": "#16a085",
              "forbidden": "0.6", "M2 (slow)": "0.6", "ground state": "k"}


def main():
    fig, ax = plt.subplots(figsize=(9.6, 8.8))

    xL, xR = 0.18, 0.74
    role_tag = {"E0": "s$\\to$E0", "M1": "s$\\to$M1", "E1": "p$\\to$E1",
                "forbidden": "p, g.s. forbidden", "M2 (slow)": "p$\\to$M2 (slow)"}
    for Ex, jp, T, w, wave, role in LEVELS:
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
            ax.text(xR + 0.175, Ex, f"[{role_tag[role]}]", va="center",
                    ha="left", fontsize=8.0, color=c, style="italic")
        else:
            ax.text(xR + 0.015, Ex, f"${jp}$  (ground state)", va="center",
                    ha="left", fontsize=10.5, color=c, fontweight="bold")

    # thresholds (labels parked at far left, staggered)
    for Eth, lab in THRESHOLDS:
        ax.axhline(Eth, xmin=0.16, xmax=0.74, color="darkorange", ls="--",
                   lw=1.2, zorder=2)
        ax.text(0.02, Eth, f"{lab} ({Eth:.2f})", va="center",
                ha="left", fontsize=8.2, color="darkorange")

    # the two teaching annotations
    ax.annotate("lowest $1^+$ is here ($\\sim$28 MeV): $\\sim$7.7 MeV ABOVE\n"
                "threshold, very broad. $\\Rightarrow$ the M1 capture is DIRECT,\n"
                "through the non-resonant $^3S_1$ continuum at\n"
                "threshold, NOT through a $1^+$ level.",
                xy=(xR + 0.05, 28.31), xytext=(0.46, 13.5),
                fontsize=8.4, color="#2471a3", ha="left",
                arrowprops=dict(arrowstyle="->", color="#2471a3"))
    ax.annotate("the ONLY sub-threshold state.\n"
                "Its tail boosts E0 capture at low $E_n$.",
                xy=(xL + 0.05, 20.21), xytext=(0.46, 8.5),
                fontsize=8.4, color="#c0392b", ha="left",
                arrowprops=dict(arrowstyle="->", color="#c0392b"))

    ax.set_ylim(-1.5, 30.5)
    ax.set_xlim(0.0, 1.18)
    ax.set_xticks([])
    ax.set_ylabel("excitation energy  $E_x$  [MeV]")
    ax.set_title("$^4$He levels near the $n+{}^3$He threshold (TUNL A=4) --- "
                 "bands = widths,\ncolour = how the state decays to the "
                 "$0^+$ ground state",
                 fontsize=11.5)

    from matplotlib.lines import Line2D
    leg = [Line2D([], [], color="#c0392b", lw=3, label="E0 $\\to$ 100% pairs ($\\gamma$-dark)"),
           Line2D([], [], color="#2471a3", lw=3, label="M1 $\\to$ $\\gamma$ + 0.35% pairs"),
           Line2D([], [], color="#16a085", lw=3, label="E1 $\\to$ $\\gamma$ + 0.35% pairs"),
           Line2D([], [], color="0.6", lw=3, label="forbidden / slow (spectator)"),
           Line2D([], [], color="darkorange", ls="--", lw=1.2, label="particle threshold")]
    ax.legend(handles=leg, loc="lower left", fontsize=8.2, framealpha=0.95)

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_he4_levels.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved -> {OUT}/fig_he4_levels.[pdf,png]")


if __name__ == "__main__":
    main()
