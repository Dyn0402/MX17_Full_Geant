#!/usr/bin/env python3
"""
make_capture_breakdown_fig.py

A schematic of what happens when a neutron is captured by 3He, drawn to make
the multipole/channel structure unambiguous: the E0 and M1 are DIFFERENT
capture channels (set by how the n and 3He spins line up), not two decay modes
of one state.  A 0+ -> 0+ transition can only be E0 (gamma forbidden).

Output: docs/e0_branch/figs/fig_capture_breakdown.{pdf,png}
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = Path(__file__).resolve().parent.parent / "docs" / "e0_branch" / "figs"


def box(ax, x, y, w, h, text, ec="0.3", fc="white", fs=10, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                 boxstyle="round,pad=0.08", ec=ec, fc=fc, lw=1.6))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", color="0.1")


def arrow(ax, x0, y0, x1, y1, color="0.4", lw=1.8, label=None, lx=0, ly=0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                 mutation_scale=16, color=color, lw=lw,
                 shrinkA=2, shrinkB=2))
    if label:
        ax.text((x0 + x1) / 2 + lx, (y0 + y1) / 2 + ly, label, fontsize=8.5,
                color=color, ha="center", style="italic")


def main():
    fig, ax = plt.subplots(figsize=(10.5, 8.2))
    ax.set_xlim(0, 12); ax.set_ylim(0, 11); ax.axis("off")

    RED, BLUE, GREY = "#c0392b", "#2471a3", "0.45"

    box(ax, 6, 10.3, 3.0, 0.8, r"$n + {}^{3}$He", fc="#eef3f8", bold=True)
    arrow(ax, 6, 9.9, 6, 9.2, label="fuse / capture", lx=1.4)
    box(ax, 6, 8.7, 4.6, 0.9,
        r"${}^{4}$He$^*$  compound  ($\approx$20.6 MeV)", fc="#eef3f8", bold=True)

    # two fates
    arrow(ax, 5.0, 8.25, 2.4, 7.05, color=GREY, lw=3.0)
    box(ax, 2.4, 6.5, 3.6, 1.2,
        "p + ${}^{3}$H\n$\\approx$100%  —  the (n,p)t\nreaction.  NO pair.",
        ec=GREY, fc="#f3f3f3")
    arrow(ax, 7.0, 8.25, 8.4, 7.15, color="0.3")
    box(ax, 8.4, 6.55, 4.6, 1.15,
        "drop to ${}^{4}$He ground state (0$^+$)\nrare, electromagnetic\n"
        "$\\sim$10$^{-8}$–10$^{-4}$ of captures", ec="0.3", fc="#fbf6ee")

    # split into the two EM channels
    arrow(ax, 7.2, 6.0, 5.3, 4.9, color=RED)
    arrow(ax, 9.6, 6.0, 10.2, 4.9, color=BLUE)

    box(ax, 4.6, 4.35, 4.2, 1.05,
        "spins OPPOSITE ($^1S_0$)\n$\\Rightarrow$ ${}^{4}$He$^*$ is 0$^+$",
        ec=RED, fc="#fdecea")
    box(ax, 10.0, 4.35, 3.4, 1.05,
        "spins ALIGNED ($^3S_1$)\n$\\Rightarrow$ ${}^{4}$He$^*$ is 1$^+$",
        ec=BLUE, fc="#eaf2f8")

    arrow(ax, 4.6, 3.8, 4.6, 3.05, color=RED,
          label="0$^+\\to$0$^+$", lx=-1.15)
    arrow(ax, 10.0, 3.8, 10.0, 3.05, color=BLUE,
          label="1$^+\\to$0$^+$", lx=1.05)

    box(ax, 4.6, 2.5, 4.2, 1.0,
        "E0 transition\n$\\gamma$ is FORBIDDEN", ec=RED, fc="#fdecea", bold=True)
    box(ax, 10.0, 2.5, 3.4, 1.0,
        "M1 transition\n$\\gamma$ allowed", ec=BLUE, fc="#eaf2f8", bold=True)

    arrow(ax, 4.6, 2.0, 4.6, 1.35, color=RED)
    arrow(ax, 10.0, 2.0, 10.0, 1.35, color=BLUE)
    box(ax, 4.6, 0.85, 4.2, 0.95,
        "100% e$^+$e$^-$ PAIR\n(the only outlet)", ec=RED, fc="#f9d5cf", bold=True)
    box(ax, 10.0, 0.85, 3.4, 0.95,
        "$\\gamma$-ray 99.65%\n+ e$^+$e$^-$ 0.35% (IPC)", ec=BLUE, fc="#d4e3f0")

    # key message banners
    ax.text(6, 5.55, "the SAME 0$\\to$0 rule that forbids the photon\n"
            "also forbids M1 here — so 0$^+$ can ONLY go by E0",
            ha="center", va="center", fontsize=9.0, color=RED, style="italic")
    ax.text(6.0, 0.1,
            "Lower $E_n$  $\\Rightarrow$  more s-wave  $\\Rightarrow$  the 0$^+$(E0) "
            "channel grows (boosted by the sub-threshold 20.21 MeV 0$^+$); "
            "the p-wave 0$^-$/2$^-$ fade out.",
            ha="center", va="bottom", fontsize=9, color="0.25")

    ax.set_title("What makes the e$^+$e$^-$ pairs: E0 and M1 are different "
                 "capture channels, not two decays of one state",
                 fontsize=11.5, fontweight="bold")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_capture_breakdown.{ext}", bbox_inches="tight",
                    dpi=150)
    plt.close(fig)
    print(f"Saved -> {OUT}/fig_capture_breakdown.[pdf,png]")


if __name__ == "__main__":
    main()
