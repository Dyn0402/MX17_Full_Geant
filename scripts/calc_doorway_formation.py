#!/usr/bin/env python3
"""
calc_doorway_formation.py

Illustrative calculation of the RELATIVE probability of FORMING the compound
4He* through its different doorway states (0+, 1+, 0-, 2-, 1-) as a function of
incoming neutron energy E_n, for n + 3He.  This is the ENTRANCE channel only --
it does NOT include the electromagnetic matrix element that decides whether the
state then radiates to the 0+ ground state (that is Part 2, and for M1/E1 it is
measured, not computed -- see doorway_states_note.md).

The entrance ("formation") cross section factorises cleanly:

    sigma_Jpi(E_n)  ~  (1/k^2)  x  g_J  x  P_l(E)  x  M_Jpi(E)
                       geometric   spin    barrier    structure

  (1/k^2)  geometric factor (pi lambdabar^2), ~1/E_cm.  Combined with the
            s-wave penetrability P_0~k it gives the textbook 1/v capture law
            (sigma ~ 1/k^2 * k = 1/v); for p-wave ~ 1/k^2 * k^3 = k^2 (the
            barrier turn-on).  Drops out of ratios at fixed l.
  g_J      = (2J+1)/[(2*1/2+1)(2*1/2+1)] = (2J+1)/4   spin-statistical weight
            -> g(1+)=3/4, g(0+)=1/4 : UNPOLARISED n+3He makes 3x more 3S1 (1+)
               than 1S0 (0+) from spin counting alone.
  P_l(E)   centrifugal penetrability: s-wave (l=0) ~ rho (1/v-law, no barrier);
            p-wave (l=1) ~ rho^3 (barrier-suppressed, turns on toward MeV).
  M_Jpi(E) "structure modulation": = 1 for a smooth direct (non-resonant)
            channel; a Breit-Wigner shape for a resonant one.  The 0+ has a
            SUB-THRESHOLD resonance (E_r=-0.37 MeV) -> its M rises toward
            threshold (boost normalised to 1 at E_n = 1 MeV); the 1+ is direct
            (M=1, lowest 1+ level is 7.7 MeV away); 0-,2-,1- are p-wave
            resonances peaking in the sub-MeV/MeV region.

IMPORTANT HONESTY NOTE (why this differs from a naive plot): the ABSOLUTE
vertical offsets between states of DIFFERENT l carry an arbitrary equal-reduced-
width assumption and are NOT physical -- only the SHAPES (energy dependence) and
the s-wave 1+/0+ comparison (same l, so P_l and 1/v cancel) are meaningful.  For
the s-wave pair the ratio is honest:
        S(1+)/S(0+) = g(1+)/g(0+) x 1/M_0+(E) = 3 / boost_0+(E),
so the spin factor of 3 (favouring 1+) is partly cancelled by the sub-threshold
0+ boost at low E_n -> the two s-wave doorways are COMPARABLE at thermal, not an
order of magnitude apart.  A previous version used an ad-hoc flat factor for the
1+ that made it spuriously dominate; that has been removed.

Parameters: TUNL A=4 (Tilley/Weller/Hale, NPA 541 (1992) 1).
Output: docs/e0_branch/figs/fig_doorway_formation.{pdf,png}
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent.parent / "docs" / "e0_branch" / "figs"

# ── constants / kinematics ────────────────────────────────────────────────────
HBARC = 197.327          # MeV fm
MU = 0.75 * 939.565      # reduced mass mc^2 [MeV] for n + 3He  (= 3/4 m_n)
A_CH = 4.0               # channel radius [fm]
S_N = 20.578             # n + 3He threshold = S_n(4He) [MeV]
EN_NORM = 1.0            # E_n [MeV] at which resonance modulations M are set = 1

# states: (label, E_x[MeV], Gamma[MeV], l, J, kind, color, linestyle)
#   kind "direct" = non-resonant direct capture (no nearby level): M=1
#   kind "sub"    = sub-threshold resonance (0+): M = boost (>1 toward threshold)
#   kind "res"    = above-threshold p-wave resonance: M peaks at E_r
STATES = [
    ("$1^+$ ($^3S_1$, s-wave, DIRECT $\\to$M1)", 28.31, 9.00, 0, 1, "direct", "#2471a3", "--"),
    ("$0^+$ (20.21, s-wave, sub-thr. $\\to$E0)", 20.21, 0.50, 0, 0, "sub",    "#c0392b", "-"),
    ("$0^-$ (21.01, p-wave $\\to$ g.s. forbidden)", 21.01, 0.84, 1, 0, "res",  "0.55",    "-"),
    ("$2^-$ (21.84, p-wave $\\to$M2)",            21.84, 2.01, 1, 2, "res",    "#16a085", "-"),
    ("$1^-$ (25.95, p-wave, T=1 $\\to$E1/GDR)",   25.95, 12.66, 1, 1, "res",   "#8e44ad", "-"),
]


def Ecm_of_En(En):
    return 0.75 * En


def k_of_Ecm(Ecm):
    return np.sqrt(2.0 * MU * np.maximum(Ecm, 1e-30)) / HBARC


def penetrability(l, rho):
    r2 = rho * rho
    if l == 0:
        return rho
    if l == 1:
        return rho * r2 / (1.0 + r2)
    if l == 2:
        return rho * r2 * r2 / (9.0 + 3.0 * r2 + r2 * r2)
    raise ValueError(l)


def bw(Ecm, Er, G):
    """Breit-Wigner shape, peak = 1 at E_r."""
    return (G / 2.0) ** 2 / ((Ecm - Er) ** 2 + (G / 2.0) ** 2)


def modulation(kind, Ecm, Er, G):
    """structure factor M(E): direct=1; sub-threshold boost normalised to 1 at
    E_n=EN_NORM; above-threshold resonance = peak-1 Breit-Wigner."""
    if kind == "direct":
        return np.ones_like(Ecm)
    if kind == "sub":
        return bw(Ecm, Er, G) / bw(Ecm_of_En(EN_NORM), Er, G)
    return bw(Ecm, Er, G)


def strength(En, Ex, G, l, J, kind):
    Ecm = Ecm_of_En(En)
    Er = Ex - S_N
    rho = k_of_Ecm(Ecm) * A_CH
    P = np.array([penetrability(l, r) for r in rho])
    gJ = (2 * J + 1) / 4.0
    # geometric factor (pi lambdabar^2) ~ 1/k^2 ~ 1/E_cm.  Combined with the
    # s-wave penetrability P_0 ~ k this gives the textbook 1/v capture law
    # (sigma ~ 1/k^2 * k = 1/k ~ 1/v); p-wave ~ 1/k^2 * k^3 ~ k^2 (barrier).
    flux = 1.0 / np.maximum(Ecm, 1e-30)
    M = modulation(kind, Ecm, Er, G)
    return flux * gJ * P * M, gJ, P, M


def main():
    En = np.logspace(-3, 1, 1400)        # 1 meV - 10 MeV
    curves, parts = {}, {}
    for (label, Ex, G, l, J, kind, c, ls) in STATES:
        S, gJ, P, M = strength(En, Ex, G, l, J, kind)
        curves[label] = S
        parts[label] = dict(gJ=gJ, P=P, M=M, l=l, kind=kind)

    # ── figure: top = doorway strength, bottom = penetrabilities ─────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.0, 8.4), sharex=True,
                                   gridspec_kw={"height_ratios": [2.15, 1],
                                                "hspace": 0.08})

    for (label, Ex, G, l, J, kind, c, ls) in STATES:
        ax1.loglog(En, curves[label], color=c, lw=2.4, ls=ls, label=label)
    ax1.axvspan(0.2, 2.0, color="gold", alpha=0.2, zorder=0)  # 0.2-2 MeV window
    ax1.text(0.6, 0.93,
             "s-wave ($0^+,1^+$): no barrier $\\Rightarrow$ rise as $1/v$, dominate low $E_n$\n"
             "p-wave ($0^-,2^-,1^-$): barrier-suppressed $\\Rightarrow$ peak sub-MeV/MeV\n"
             "$0^+$ sub-threshold boost makes it COMPARABLE to $1^+$ at thermal\n"
             "(the $3{:}1$ spin weight favouring $1^+$ is largely cancelled)",
             transform=ax1.transAxes, va="top", ha="left", fontsize=7.8,
             bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.92))
    ax1.text(0.015, 0.03,
             "ABSOLUTE heights across different $\\ell$ assume equal reduced "
             "widths (arbitrary);\nonly SHAPES + the same-$\\ell$ $1^+/0^+$ "
             "comparison are physical. EM (M1/E1/E0)\nstrength is a SEPARATE "
             "factor, not shown (measured for M1/E1 -- see text).",
             transform=ax1.transAxes, va="bottom", ha="left", fontsize=7.2,
             color="0.35", style="italic")
    ax1.set_ylabel("formation cross section  [arb., $\\propto\\sigma$]")
    ax1.set_title("Which $^4$He$^*$ doorway do we FORM, vs neutron energy?\n"
                  "$\\sigma \\propto (1/k^2)\\times g_J\\times P_\\ell(E)\\times M(E)$  "
                  "(s-wave $\\to 1/v$; entrance channel only; TUNL A=4, $a=4$ fm)")
    ax1.legend(loc="lower left", fontsize=8.0, framealpha=0.95)
    ax1.grid(alpha=0.3, which="both")
    ymax = max(np.nanmax(v) for v in curves.values())
    ax1.set_ylim(ymax * 3e-5, ymax * 5)
    ax1.set_xlim(1e-3, 1e1)

    # penetrabilities
    for l, c, lab in [(0, "#2471a3", r"$P_0$ (s-wave: $0^+,1^+$)"),
                      (1, "C1", r"$P_1$ (p-wave: $0^-,2^-,1^-$)"),
                      (2, "C2", r"$P_2$ (d-wave)")]:
        rho = k_of_Ecm(Ecm_of_En(En)) * A_CH
        P = np.array([penetrability(l, r) for r in rho])
        ax2.loglog(En, P, color=c, lw=1.8, label=lab)
    ax2.set_xlabel("incoming neutron energy  $E_n$  [MeV]")
    ax2.set_ylabel("penetrability $P_\\ell$")
    ax2.set_ylim(1e-6, 5)
    ax2.legend(loc="lower right", fontsize=9)
    ax2.grid(alpha=0.3, which="both")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_doorway_formation.{ext}", bbox_inches="tight",
                    dpi=150)
    plt.close(fig)

    # ── console: the 1+/0+ decomposition (the user's question) ───────────────
    lab1p, lab0p, lab0m = STATES[0][0], STATES[1][0], STATES[2][0]
    s1p, s0p, s0m = curves[lab1p], curves[lab0p], curves[lab0m]
    print("Decomposing S(1+)/S(0+)  [both s-wave: 1/v and P_0 cancel exactly]")
    print(f"  spin weight ratio g(1+)/g(0+) = {parts[lab1p]['gJ']/parts[lab0p]['gJ']:.1f}  (FIXED, physical)")
    print(f"  {'E_n[MeV]':>9} {'boost_0+':>9} {'S(1+)/S(0+)':>12} {'S(0+)/S(0-)':>12}")
    for En_pt in [1e-3, 1e-2, 0.1, 0.3, 0.6, 1.0, 3.0]:
        j = np.argmin(np.abs(En - En_pt))
        boost = parts[lab0p]['M'][j]
        print(f"  {En_pt:>9.3g} {boost:>9.2f} {s1p[j]/s0p[j]:>12.2f} {s0p[j]/s0m[j]:>12.2e}")
    i_pk = np.argmax(s0m)
    print(f"\n0- (p-wave) doorway peaks at E_n = {En[i_pk]*1e3:.0f} keV "
          f"(E_cm = {Ecm_of_En(En[i_pk])*1e3:.0f} keV)")
    # s/p crossover for the 0+ (s) vs 0- (p)
    sign = np.sign(s0p / s0m - 1.0)
    xover = En[np.where(np.diff(sign) != 0)[0]]
    if len(xover):
        print(f"s-wave 0+ vs p-wave 0- crossover at E_n ~ {xover[-1]*1e3:.0f} keV")
    print("\nTakeaway: the 1+ does NOT dominate by an order of magnitude. At "
          "thermal\n  S(1+)/S(0+) ~ 3/boost ~ O(1); the gap only opens toward "
          "MeV as the 0+\n  sub-threshold boost fades. The radiative M1 vs E0 "
          "hierarchy is separate.")
    print(f"\nSaved -> {OUT}/fig_doorway_formation.[pdf,png]")


if __name__ == "__main__":
    main()
