#!/usr/bin/env python3
"""
calc_direct_vs_resonant.py

How much of each capture channel is "resonant" vs "direct" (continuum)?

A level acts as a sharp resonant doorway only if its width Gamma is SMALLER than
its distance |E_r| from the capture window.  For 4He above the n+3He threshold
the levels are very broad (Gamma ~ a few MeV), comparable to or larger than
their distance from threshold -- so most of the capture is effectively DIRECT
(continuum) capture, modulated by broad resonant structures, not sharp
resonant capture.

We quantify this two ways:
  (1) Gamma/|E_r| for each level (>~1 => not a sharp resonance => ~direct).
  (2) the Breit-Wigner "resonant factor" g(E) = (Gamma/2)^2/((E-E_r)^2+(Gamma/2)^2)
      across the capture window: flat => direct background; peaked => resonant.

Key results printed: the subthreshold 0+ gives a real (but broad) low-energy
ENHANCEMENT of the E0 capture (~10x at threshold vs 1 MeV), while the 1+ (M1)
is essentially flat = pure direct continuum capture.

Output: docs/e0_branch/figs/fig_direct_vs_resonant.{pdf,png}
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent.parent / "docs" / "e0_branch" / "figs"
S_N = 20.578   # n+3He threshold [MeV]

# (label, E_x, Gamma, color, g.s. transition)
STATES = [
    ("$1^+$ ($^3S_1$, M1)", 28.31, 9.00, "#2471a3", "M1"),
    ("$0^+$ (20.21, E0)",   20.21, 0.50, "#c0392b", "E0"),
    ("$0^-$ (21.01)",       21.01, 0.84, "C1",      "—"),
    ("$2^-$ (21.84)",       21.84, 2.01, "#16a085", "M2"),
    ("$1^-$ (24.25, E1)",   24.25, 6.10, "#8e44ad", "E1"),
]


def g_factor(Ecm, Er, G):
    return (G / 2.0) ** 2 / ((Ecm - Er) ** 2 + (G / 2.0) ** 2)


def main():
    En = np.logspace(-3, 1.0, 1500)     # MeV
    Ecm = 0.75 * En

    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    for label, Ex, G, c, gs in STATES:
        Er = Ex - S_N
        g = g_factor(Ecm, Er, G)
        ls = "--" if "1^+" in label else "-"
        ax.loglog(En, g, color=c, lw=2.3, ls=ls,
                  label=f"{label}  $\\to${gs},  $\\Gamma/|E_r|$={G/abs(Er):.1f}")
    ax.axvspan(0.2, 2.0, color="gold", alpha=0.2, zorder=0)
    ax.text(0.63, 1.25, "0.2–2 MeV\nwindow", ha="center", fontsize=8.5,
            color="0.4")
    ax.axhline(1.0, color="0.6", ls=":", lw=1)
    ax.text(1.3e-3, 1.06, "on resonance", fontsize=7.5, color="0.5")
    ax.text(1.3e-3, 0.018, "flat $\\Rightarrow$ effectively DIRECT (continuum) capture;\n"
            "peaked $\\Rightarrow$ resonant", fontsize=8.3, color="0.25",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    ax.set_xlabel("incoming neutron energy  $E_n$  [MeV]")
    ax.set_ylabel("Breit–Wigner resonant factor  $g(E)$  (peak $=1$)")
    ax.set_title("Resonant vs direct character of each $^4$He capture doorway\n"
                 "(4He levels are broad: $\\Gamma\\gtrsim|E_r|$ for all "
                 "$\\Rightarrow$ capture is mostly direct)")
    ax.set_ylim(1e-2, 2.0)
    ax.set_xlim(1e-3, 1e1)
    ax.legend(loc="upper left", fontsize=8.2, ncol=1, framealpha=0.95)
    ax.grid(alpha=0.3, which="both")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_direct_vs_resonant.{ext}", bbox_inches="tight",
                    dpi=150)
    plt.close(fig)

    # ── numbers ──────────────────────────────────────────────────────────────
    print(f"{'state':>22} {'E_r(cm)':>9} {'Gamma':>6} {'G/|Er|':>7} "
          f"{'g(thr)':>7} {'g(1MeVcm)':>9} {'thr/1MeV':>9}  character")
    for label, Ex, G, c, gs in STATES:
        Er = Ex - S_N
        gthr = g_factor(0.0, Er, G)
        g1 = g_factor(1.0, Er, G)
        ratio = gthr / g1
        if G / abs(Er) > 3 or abs(Er) > 3:
            char = "direct (flat)"
        elif Er < 0:
            char = "subthreshold boost"
        else:
            char = f"resonant, peak En~{Er/0.75*1e3:.0f}keV"
        lab = label.split("(")[0].strip()
        print(f"{lab:>22} {Er:>9.2f} {G:>6.2f} {G/abs(Er):>7.1f} "
              f"{gthr:>7.3f} {g1:>9.3f} {ratio:>9.2f}  {char}")
    print("\nReadout:")
    print("  1+ (M1): g nearly flat (thr/1MeV ~ 0.8) -> PURE DIRECT continuum "
          "capture; the 54 ub is direct.")
    print("  0+ (E0): subthreshold 0+ enhances capture ~10x at threshold vs "
          "1 MeV -> a real (broad) resonant boost at low E_n.")
    print("  0-,2-,1-: broad p-wave resonances; matter only toward MeV.")
    print(f"\nSaved -> {OUT}/fig_direct_vs_resonant.[pdf,png]")


if __name__ == "__main__":
    main()
