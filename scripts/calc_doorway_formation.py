#!/usr/bin/env python3
"""
calc_doorway_formation.py

Illustrative calculation of the RELATIVE probability of forming the compound
4He* through the 0+ (20.21 MeV) vs 0- (21.01 MeV) doorway states as a function
of incoming neutron energy E_n, for n + 3He.

The physics in one line: a 0+ state (parity +) is reached by s-wave (l=0)
neutron capture, a 0- state (parity -) by p-wave (l=1).  s-wave has no
centrifugal barrier -> turns on as 1/v at threshold; p-wave is barrier-
suppressed at low energy -> vanishes ~v at threshold and rises with E_n.  On
top of that sit the resonance energies (0+ is 0.37 MeV BELOW the n+3He
threshold = subthreshold; 0- is 0.43 MeV above).

Single-level Breit-Wigner with neutral-particle penetrabilities:

    S_J(E) ∝ (1/k^2) g_J  Γ_n^J(E)  Γ_J / [ (E_cm - E_r)^2 + (Γ_J/2)^2 ]
    Γ_n^J(E) = 2 P_l(E) γ^2 ,   γ^2 = θ^2 (ħc)^2 / (μc^2 a^2)

with equal dimensionless reduced widths θ^2 for the comparison.  This is an
illustration of the ENTRANCE-channel / doorway strength -- it shows WHY 0+
favours low E_n and 0- needs ~MeV.  It is NOT the absolute rate or the EM/pair
branching (those need the outgoing widths and the real multi-level 4He
R-matrix, Hale).  Parameters: TUNL A=4 (Tilley/Weller/Hale, NPA 541 (1992) 1).

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
S_N = 20.578            # n + 3He threshold = S_n(4He) [MeV]
THETA2 = 0.1            # dimensionless reduced width (Wigner units)

# states: (label, E_x[MeV], Gamma[MeV], l, J, kind)
#   kind "res"    = resonant doorway (Breit-Wigner around E_r)
#   kind "direct" = non-resonant direct capture (no nearby level): the 1+ M1
#                   capture goes through the 3S1 continuum, since the lowest 1+
#                   level is ~28 MeV away.  Modelled as a reference-width
#                   "resonance" so it is just the smooth s-wave 1/v entrance.
G_REF = 1.0   # MeV, reference width for the direct (non-resonant) 1+ channel
STATES = [
    ("$1^+$ ($^3S_1$, s, direct)",        28.31, G_REF, 0, 1, "direct"),
    ("$0^+$ (20.21, s, sub-thr.)",        20.21, 0.50,  0, 0, "res"),
    ("$0^-$ (21.01, p-wave)",             21.01, 0.84,  1, 0, "res"),
    ("$2^-$ (21.84, p-wave)",             21.84, 2.01,  1, 2, "res"),
]
COLORS = ["#2471a3", "#c0392b", "C1", "#16a085"]


def Ecm_of_En(En):                       # CM energy from lab neutron energy
    return 0.75 * En


def k_of_Ecm(Ecm):                       # wavenumber [1/fm]
    return np.sqrt(2.0 * MU * np.maximum(Ecm, 1e-30)) / HBARC


def penetrability(l, rho):
    """Neutral-particle penetrability P_l(rho), rho = k a."""
    r2 = rho * rho
    if l == 0:
        return rho
    if l == 1:
        return rho * r2 / (1.0 + r2)
    if l == 2:
        return rho * r2 * r2 / (9.0 + 3.0 * r2 + r2 * r2)
    raise ValueError(l)


def main():
    En = np.logspace(-3, 1, 1200)        # 1 meV - 10 MeV
    Ecm = Ecm_of_En(En)
    k = k_of_Ecm(Ecm)
    rho = k * A_CH
    gamma2 = THETA2 * HBARC**2 / (MU * A_CH**2)     # reduced width [MeV]

    curves = {}
    for (label, Ex, G, l, J, kind) in STATES:
        Er = Ex - S_N                                # resonance E in CM
        P = np.array([penetrability(l, r) for r in rho])
        gJ = (2 * J + 1) / 4.0                        # (2J+1)/[(2*1/2+1)^2]
        Gn = 2.0 * P * gamma2                         # energy-dependent neutron width
        if kind == "direct":
            # non-resonant: flat reference factor -> just the s-wave 1/v entrance
            F = G / ((G / 2.0) ** 2)                  # = 4/G_ref, energy-independent
        else:
            F = G / ((Ecm - Er) ** 2 + (G / 2.0) ** 2)
        S = (1.0 / k**2) * gJ * Gn * F
        curves[label] = S

    # ── figure: top = doorway strength, bottom = penetrabilities ─────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.6, 8.0), sharex=True,
                                   gridspec_kw={"height_ratios": [2, 1],
                                                "hspace": 0.08})

    for (label, *_), c in zip(STATES, COLORS):
        ls = "--" if label.startswith("$1^+$") else "-"
        ax1.loglog(En, curves[label], color=c, lw=2.4, ls=ls, label=label)
    ax1.axvspan(2e5 * 1e-6, 2e6 * 1e-6, color="gold", alpha=0.2, zorder=0)  # 0.2-2 MeV
    ax1.text(0.06, 0.94, "s-wave channels ($0^+,1^+$): rise as $1/v$, dominate "
             "low $E_n$\np-wave ($0^-,2^-$): barrier-suppressed, peak in "
             "sub-MeV/MeV", transform=ax1.transAxes, va="top", fontsize=8.2,
             bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    ax1.set_ylabel("relative formation strength  [arb.]")
    ax1.set_title("Which $^4$He$^*$ doorway do we form, vs neutron energy?\n"
                  "(single-level Breit–Wigner $\\times$ penetrability; TUNL A=4, "
                  "$a=4$ fm. $1^+$ is direct/non-resonant.)")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(alpha=0.3, which="both")
    ymax = max(np.nanmax(v) for v in curves.values())
    ax1.set_ylim(ymax * 1e-5, ymax * 4)

    # penetrabilities
    for l, c, lab in [(0, "#2471a3", r"$P_0$ (s-wave: $0^+,1^+$)"),
                      (1, "C1", r"$P_1$ (p-wave: $0^-,2^-$)"),
                      (2, "C2", r"$P_2$ (d-wave)")]:
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

    # ── console: ratios + crossover ──────────────────────────────────────────
    s1p = curves[STATES[0][0]]   # 1+ (direct s-wave)
    s0p = curves[STATES[1][0]]   # 0+ (sub-threshold s-wave)
    s0m = curves[STATES[2][0]]   # 0- (p-wave)
    ratio = s0p / s0m            # s-wave 0+ vs p-wave 0-
    i_pk = np.argmax(s0m)
    sign = np.sign(ratio - 1.0)
    xover = En[np.where(np.diff(sign) != 0)[0]]
    print(f"reduced width gamma^2 = {gamma2:.3f} MeV (theta^2={THETA2}, a={A_CH} fm)")
    print(f"0- doorway peaks at E_n = {En[i_pk]*1e3:.0f} keV "
          f"(E_cm = {Ecm_of_En(En[i_pk])*1e3:.0f} keV)")
    print(f"{'E_n':>9} {'S(0+)/S(0-)':>12} {'S(0+)/S(1+)':>12}  (s/p, and s-wave pair)")
    for En_pt in [1e-3, 1e-2, 0.1, 0.3, 0.6, 1.0, 3.0]:
        j = np.argmin(np.abs(En - En_pt))
        print(f"  {En_pt:>7.3g} {ratio[j]:>12.2e} {s0p[j]/s1p[j]:>12.2e}")
    if len(xover):
        print(f"s-wave/p-wave crossover S(0+)=S(0-) at E_n ~ {xover[0]*1e3:.0f} keV")
    print("note: S(0+)/S(1+) is ~flat (both s-wave); its absolute value is the\n"
          "      unknown E0/M1 matrix-element ratio -- shapes robust, height not.")
    print(f"\nSaved -> {OUT}/fig_doorway_formation.[pdf,png]")


if __name__ == "__main__":
    main()
