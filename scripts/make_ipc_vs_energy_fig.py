#!/usr/bin/env python3
"""
make_ipc_vs_energy_fig.py

Phase-0 deliverable for the IPC roadmap: the number of internal-pair (e+e-)
pairs produced per day as a function of neutron energy, decomposed into the
piece we can compute now (term 1: the M1/E1 photon-emitting radiative capture,
times the high-energy IPC coefficient) and the piece that is still missing
(term 2: the gamma-dark E0 channel, drawn as a labelled gap).

    N_pairs(E_n)/pulse = N_beam(E_n) * [ P_radcap(E_n)*alpha_IPC   (term 1, HAVE)
                                       + P_E0(E_n)*1               (term 2, MISSING) ]

Term 1 uses the directly-counted 3He(n,gamma)4He rate per decade from the
5e8-neutron campaign (analysis/mev/mev_rates.json) and the high-energy IPC
coefficient alpha_IPC.  alpha_IPC at ~20.6 MeV is NOT in BrIcc (its tables stop
at 6 MeV); it is anchored to measured high-energy transitions of the same
energy/multipolarity:
    12C 15.1 MeV M1 : alpha_pi = 3.3(5)e-3  (measured)
    8Be 18.15 MeV M1: ~3-4e-3
=> alpha_IPC(20.6 MeV) ~ 3.5e-3, band 2.5-4.5e-3.  (Alberto's table used
2.1e-3 -- shown for comparison.)

Output: docs/e0_branch/figs/fig_ipc_vs_energy.{pdf,png}
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RATES = ROOT / "analysis" / "mev" / "mev_rates.json"
OUT = ROOT / "docs" / "e0_branch" / "figs"

ALPHA_ANCHOR = 3.5e-3      # high-energy IPC coefficient, 12C/8Be-anchored
ALPHA_LO, ALPHA_HI = 2.5e-3, 4.5e-3
ALPHA_TABLE = 2.1e-3       # Alberto's rate-table value, for comparison


def main():
    d = json.load(open(RATES))
    ppd = d["pulses_per_day"]
    rows = d["decades"]
    Elo = np.array([r["E_lo_eV"] for r in rows])
    Ehi = np.array([r["E_hi_eV"] for r in rows])
    Ecen = np.sqrt(Elo * Ehi)
    rad_pp = np.array([r["rad_per_pulse"] for r in rows])        # (n,g)/pulse/decade
    n_rad = np.array([r["n_rad"] for r in rows])
    lo68 = np.array([r["n_rad_68"][0] for r in rows])
    hi68 = np.array([r["n_rad_68"][1] for r in rows])
    # scale Poisson band on counts to the rate
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(n_rad > 0, rad_pp / n_rad, rad_pp / 1.0)
    # for zero-count decades rad_pp=0; use the 68% upper as an upper limit via weight
    w = d["weight_per_neutron"]
    ipc_day = rad_pp * ALPHA_ANCHOR * ppd
    ipc_day_lo = lo68 * w * ALPHA_LO * ppd
    ipc_day_hi = hi68 * w * ALPHA_HI * ppd
    ipc_day_table = rad_pp * ALPHA_TABLE * ppd

    fig, ax = plt.subplots(figsize=(9, 5.6))

    # band from alpha + Poisson
    ax.fill_between(Ecen, np.maximum(ipc_day_lo, 1e-6), ipc_day_hi,
                    color="C0", alpha=0.2, step="mid",
                    label=r"term 1 band ($\alpha_{\rm IPC}=2.5$–$4.5\times10^{-3}$, "
                          "Poisson)")
    has = n_rad > 0
    ax.plot(Ecen[has], ipc_day[has], "o-", color="C0", ms=7, lw=2,
            label=r"term 1: M1/E1 IPC ($\alpha_{\rm IPC}=3.5\times10^{-3}$)")
    # upper limits where zero counts
    for i in np.where(~has)[0]:
        ax.annotate("", xy=(Ecen[i], ipc_day_hi[i]),
                    xytext=(Ecen[i], ipc_day_hi[i] * 3),
                    arrowprops=dict(arrowstyle="-|>", color="C0", alpha=0.6))
    ax.plot(Ecen[has], ipc_day_table[has], "s--", color="0.5", ms=5, lw=1.2,
            label=r"term 1 with table $\alpha=2.1\times10^{-3}$")

    # E0 term: drawn as an unknown gap riding on top of term 1
    ax.fill_between(Ecen[has], ipc_day[has], ipc_day[has] * 2.0,
                    color="C3", alpha=0.18, hatch="//", step=None,
                    label="term 2: E0 (γ-dark) — magnitude TBD (nuclear input)")
    ax.annotate("E0 channel:\nabsent from ENDF + Geant4 + generator\n"
                "size set by $\\rho^2$(E0) — needs nuclear input",
                xy=(3e6, ipc_day[has][-2] * 1.5 if has.any() else 1),
                xytext=(2e3, 30), fontsize=8.5, color="C3",
                arrowprops=dict(arrowstyle="->", color="C3", alpha=0.7))

    ax.axvspan(2e5, 2e6, color="gold", alpha=0.22, zorder=0)
    ax.text(8.5e5, ax.get_ylim()[1] if False else 200, "0.2–2 MeV\nwindow",
            ha="center", fontsize=8.5, color="0.4")

    tot = (rad_pp * ALPHA_ANCHOR * ppd).sum()
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"neutron energy  $E_n$  [eV]")
    ax.set_ylabel("IPC pairs produced / day / decade")
    ax.set_title("Internal-pair (e$^+$e$^-$) production vs neutron energy\n"
                 f"term 1 total $\\approx${tot:.0f} pairs/day "
                 r"($5\times10^{8}$-neutron campaign); E0 term still to be added")
    ax.set_ylim(1e-2, 1e3)
    ax.set_xlim(1e-1, 1e8)
    ax.legend(fontsize=8.5, loc="upper left")
    ax.grid(alpha=0.3, which="both")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_ipc_vs_energy.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"term-1 total (alpha={ALPHA_ANCHOR}): {tot:.1f} IPC pairs/day")
    print(f"  (table alpha={ALPHA_TABLE}: {(rad_pp*ALPHA_TABLE*ppd).sum():.1f}/day)")
    print(f"Saved -> {OUT}/fig_ipc_vs_energy.[pdf,png]")


if __name__ == "__main__":
    main()
