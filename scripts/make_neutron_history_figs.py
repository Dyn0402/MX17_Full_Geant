#!/usr/bin/env python3
"""
make_neutron_history_figs.py

Supporting figures for docs/neutron_histories/neutron_histories_note.md
(Alberto's question 1: are neutron scattering / slowing-down / capture in the
gas cell included, and could in-gas moderation enhance the excitation rate?).

Inputs (local, from the lxplus MeV campaign scan):
    analysis/mev/mev_rates.json     -- per-decade (n,p) and direct (n,g) budget

Outputs:
    docs/neutron_histories/figs/fig_excitation_fraction.pdf
    docs/neutron_histories/figs/fig_excitation_fraction.png

The figure compares the Geant4 full-transport 3He(n,p)t absorption fraction
per decade (= number of 4He* excitations per incident beam neutron) against a
naive *single-pass on-axis* thin-target estimate, and annotates the elastic
moderation scale.  It is NOT a claim of enhancement: the point is that the
full transport sits at/below the on-axis single pass (beam wider than the
2 cm bore => beam-averaged chord < on-axis), and that the efficient-but-
geometrically-limited in-gas moderation is fully contained in the simulation.
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RATES = ROOT / "analysis" / "mev" / "mev_rates.json"
OUT = ROOT / "docs" / "neutron_histories" / "figs"

# ── He-3 target physics constants (from DetectorConstruction.cc) ──────────────
RHO_HE3 = 62.7e-3        # g/cm3 at 500 bar
M_HE3 = 3.0160293        # g/mol
NA = 6.02214076e23
N_HE3 = RHO_HE3 / M_HE3 * NA          # number density [1/cm3]
L_AXIS = 6.0             # on-axis gas path: cyl L=40mm + 2x hemisphere r=10mm [cm]
NL_AXIS = N_HE3 * L_AXIS * 1e-24      # areal density [atoms/barn]
SIG_NP_TH = 5333.0       # 3He(n,p)t cross section at 0.0253 eV [barn]
SIG_EL = 3.3             # representative 3He elastic cross section [barn]


def sig_np(E_eV):
    """1/v extrapolation of the (n,p)t cross section (good below ~100 keV)."""
    return SIG_NP_TH * np.sqrt(0.0253 / E_eV)


def main():
    decades = json.load(open(RATES))["decades"]
    Elo = np.array([r["E_lo_eV"] for r in decades])
    Ehi = np.array([r["E_hi_eV"] for r in decades])
    Ecen = np.sqrt(Elo * Ehi)
    g4_frac = np.array([r.get("np_frac_of_beam", 0.0) for r in decades])

    # single-pass on-axis absorption probability, 1 - exp(-nL sigma_np)
    tau_axis = NL_AXIS * sig_np(Ecen)
    P_axis = 1.0 - np.exp(-tau_axis)

    lam_el = 1.0 / (N_HE3 * SIG_EL * 1e-24)        # elastic MFP [cm]
    p_scat = 1.0 - np.exp(-L_AXIS / lam_el)        # P(>=1 elastic scatter on axis)

    fig, ax = plt.subplots(figsize=(8.6, 5.2))

    ax.step(np.r_[Elo, Ehi[-1]], np.r_[g4_frac, g4_frac[-1]], where="post",
            color="C0", lw=2.4,
            label=r"Geant4 full transport: $^3$He(n,p)t / incident neutron")
    ax.fill_between(np.r_[Elo, Ehi[-1]], 0, np.r_[g4_frac, g4_frac[-1]],
                    step="post", alpha=0.12, color="C0")
    ax.plot(Ecen, P_axis, "s--", color="C3", ms=6,
            label=r"single-pass, on-axis ($1-e^{-nL\sigma_{np}}$, $1/v$)")

    # opacity-regime shading
    ax.axvspan(1e-3, 1e3, color="0.85", alpha=0.5, zorder=0)
    ax.axvspan(2e5, 2e6, color="gold", alpha=0.25, zorder=0)
    ax.text(3e-2, 0.93, "opaque gas\n(thermal–eV)", ha="center", fontsize=8.5,
            color="0.35")
    ax.text(8.5e5, 0.93, "0.2–2 MeV\nanalysis window", ha="center",
            fontsize=8.5, color="0.45")

    txt = (r"$^3$He elastic moderation scale:" "\n"
           rf"  $\lambda_{{\rm el}}\approx{lam_el:.0f}$ cm  vs  gas path $\approx{L_AXIS:.0f}$ cm"
           "\n"
           rf"  $\Rightarrow$ P(scatter on axis) $\approx{p_scat:.2f}$" "\n"
           r"  up to 75% energy loss per collision" "\n"
           r"  ($\xi=0.54$, $\sim$2 collisions / decade)")
    ax.text(0.985, 0.04, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.3, family="monospace",
            bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9))

    ax.set_xscale("log")
    ax.set_xlabel(r"incident neutron energy  $E_n$  [eV]")
    ax.set_ylabel(r"fraction forming $^4$He$^*$  (per incident neutron)")
    ax.set_title("In-gas neutron capture (= $^4$He$^*$ excitations): "
                 "full transport vs single-pass estimate")
    ax.set_ylim(0, 1.0)
    ax.set_xlim(1e-3, 1e8)
    ax.legend(loc="center left", fontsize=9.5, framealpha=0.95)
    ax.grid(alpha=0.3, which="both")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_excitation_fraction.{ext}", bbox_inches="tight",
                    dpi=150)
    plt.close(fig)

    # console summary used in the note
    print(f"n(He3)        = {N_HE3:.3e} /cm3")
    print(f"nL_axis       = {NL_AXIS:.4f} at/b")
    print(f"elastic MFP   = {lam_el:.1f} cm   P(scatter on axis) = {p_scat:.2f}")
    print(f"{'decade':>14} {'G4 np/beam':>11} {'P_axis':>8} {'G4/P_axis':>10}")
    for lo, hi, g4, P in zip(Elo, Ehi, g4_frac, P_axis):
        print(f"{lo:.0e}-{hi:.0e} {g4:>11.3f} {P:>8.3f} "
              f"{g4 / P if P > 0 else 0:>10.2f}")
    print(f"\nSaved -> {OUT}/fig_excitation_fraction.[pdf,png]")


if __name__ == "__main__":
    main()
