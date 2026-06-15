#!/usr/bin/env python3
"""
make_e0_luminosity_fig.py

The "from cross section to statistics" bridge figure.  The X17 (or E0) yield is

    N_X / day  =  phi(E)  x  f_breakup(E)  x  [ sigma_X(E) / sigma_np(E) ]  x  pulses/day
                  flux       fraction of      X17 per breakup (cross-section
                  [/pulse]   beam that         ratio)
                             breaks up

equivalently  N_X = N_(n,p) x sigma_X/sigma_np  -- the (n,p)t breakup, directly
counted in Geant4, is the luminosity monitor (the opacity cancels in the ratio
because (n,p) is the dominant absorber, so it attenuates both channels equally).

This figure shows the four simulation-/data-derived factors, vs neutron energy:
  (a) total beam flux phi(E)            [n/pulse/decade]   -- EAR2 spectrum (sim)
  (b) fraction of beam undergoing (n,p)t in the gas        -- opaque (sub-keV)
                                                              -> thin (MeV) (sim)
  (c) effective luminosity L = N_(n,p)/sigma_np [1/barn/pulse] -- multiply
                                                              sigma_X by this
  (d) breakups per X17 = sigma_np / sigma_X17              -- the background-per-
                                                              signal burden

Inputs: analysis/mev/mev_rates.json (per-decade beam flux, np fraction, direct
(n,p)/(n,gamma) counts) and data/He3.h5 (ENDF sigma_np, sigma_ng).

Output: docs/e0_branch/figs/fig_e0_luminosity.{pdf,png}
        (copied to docs/report/figs/fig_e0_lumi.{pdf,png})
"""
from pathlib import Path
import json
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
H5 = ROOT / "data" / "He3.h5"
RATES = ROOT / "analysis" / "mev" / "mev_rates.json"
OUT = ROOT / "docs" / "e0_branch" / "figs"
REPORT = ROOT / "docs" / "report" / "figs"
TEMP = "294K"

ER_0P, G_0P = 20.21 - 20.578, 0.50
ALPHA_IPC, BR_X17 = 3.5e-3, 0.025
F_CENTRAL = 3e-3


def read_xs(rx):
    with h5py.File(H5, "r") as f:
        g = f["He3"]
        E = g[f"energy/{TEMP}"][...]
        ds = g[f"reactions/{rx}/{TEMP}/xs"]
        i0 = ds.attrs.get("threshold_idx", 0)
        xs = ds[...]
        return E[i0:i0 + xs.shape[0]], xs


def loginterp(x, xp, fp):
    return 10 ** np.interp(np.log10(x), np.log10(xp),
                           np.log10(np.clip(fp, 1e-40, None)))


def g0plus(E_cm_MeV):
    num = ER_0P**2 + (G_0P / 2)**2
    return num / ((E_cm_MeV - ER_0P)**2 + (G_0P / 2)**2)


def main():
    d = json.load(open(RATES))
    rows = d["decades"]
    Elo = np.array([r["E_lo_eV"] for r in rows])
    Ehi = np.array([r["E_hi_eV"] for r in rows])
    Ecen = np.sqrt(Elo * Ehi)
    beam = np.array([r["beam_per_pulse"] for r in rows])
    npf = np.array([r["np_frac_of_beam"] for r in rows])
    Nnp = beam * npf                                     # (n,p)/pulse/decade

    Eng, sng = read_xs("reaction_102")
    Enp, snp = read_xs("reaction_103")
    snp_dec = loginterp(Ecen, Enp, snp)                  # barn at decade centre
    L = Nnp / snp_dec                                    # 1/barn/pulse

    # continuous grid for the cross-section-ratio panel
    Eg = np.logspace(-3, 7, 1500)
    s_ng = loginterp(Eg, Eng, sng)
    s_np = loginterp(Eg, Enp, snp)
    boost = g0plus(0.75 * Eg * 1e-6)
    sx_m1e1 = s_ng * ALPHA_IPC * BR_X17
    sx_tot = sx_m1e1 + F_CENTRAL * s_ng * boost * BR_X17   # +E0 (central f)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2))
    (aA, aB), (aC, aD) = axes
    for ax in axes.ravel():
        ax.set_xscale("log")
        ax.axvspan(2e5, 2e6, color="gold", alpha=0.18, zorder=0)
        ax.axvspan(1e-3, 1e3, color="C0", alpha=0.06, zorder=0)
        ax.grid(alpha=0.3, which="both")
        ax.set_xlim(1e-3, 1e7)

    def stepc(ax, x, y, **kw):
        ax.plot(x, y, "o-", ms=5, **kw)

    # (a) flux
    aA.set_yscale("log")
    stepc(aA, Ecen, beam, color="#34495e")
    aA.set_ylabel("beam flux  $\\phi$  [n/pulse/decade]")
    aA.set_title("(a) total neutron flux (EAR2 spectrum, sim)")

    # (b) breakup fraction
    stepc(aB, Ecen, npf, color="#c0392b")
    aB.set_yscale("log")
    aB.set_ylabel("fraction of beam $\\to$ (n,p)t in gas")
    aB.set_title("(b) breakup (interaction) fraction: opaque $\\to$ thin")
    aB.text(3e-2, 0.4, "opaque\n(~70%)", fontsize=8, color="#c0392b", ha="center")
    aB.text(3e6, 0.012, "thin\n(~4%)", fontsize=8, color="#c0392b", ha="center")

    # (c) effective luminosity
    aC.set_yscale("log")
    stepc(aC, Ecen, L, color="#16a085")
    aC.set_xlabel("neutron energy  $E_n$  [eV]")
    aC.set_ylabel("$L = N_{(n,p)}/\\sigma_{np}$  [barn$^{-1}$/pulse]")
    aC.set_title("(c) effective luminosity = breakups / $\\sigma_{np}$")
    aC.text(0.03, 0.05, "multiply $\\sigma_X$[barn] by this\n(× pulses/day) "
            "$\\to$ X/day", transform=aC.transAxes, fontsize=8, color="#16a085",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    # (d) breakups per X17
    aD.set_yscale("log")
    aD.loglog(Eg, s_np / sx_m1e1, color="#2471a3", lw=2.2,
              label=r"per M1+E1 X17")
    aD.loglog(Eg, s_np / sx_tot, color="#c0392b", lw=1.8, ls="--",
              label=r"per total X17 (+E0, $f=3\times10^{-3}$)")
    aD.set_xlabel("neutron energy  $E_n$  [eV]")
    aD.set_ylabel("$\\sigma_{np}/\\sigma_{\\rm X17}$  (breakups per X17)")
    aD.set_title("(d) background burden: (n,p)t breakups per X17 produced")
    aD.legend(fontsize=8, loc="upper right")

    fig.suptitle("From cross section to statistics: "
                 r"$N_X/\mathrm{day}=\phi\times f_{\rm breakup}\times"
                 r"(\sigma_X/\sigma_{np})\times$pulses/day"
                 "   (gold = 0.2–2 MeV window, blue = sub-keV)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    OUT.mkdir(parents=True, exist_ok=True)
    REPORT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_e0_luminosity.{ext}", bbox_inches="tight", dpi=150)
        fig.savefig(REPORT / f"fig_e0_lumi.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

    # console
    print(f"{'E_n[eV]':>10} {'flux/pulse':>11} {'np_frac':>8} {'L[1/b/pulse]':>13} "
          f"{'np/X17(M1E1)':>13}")
    for i in [2, 6, 7, 8]:
        sxi = loginterp(Ecen[i], Eg, sx_m1e1)
        npx = loginterp(Ecen[i], Eg, s_np) / sxi
        print(f"{Ecen[i]:>10.2e} {beam[i]:>11.2e} {npf[i]:>8.3f} {L[i]:>13.3e} "
              f"{npx:>13.2e}")
    print(f"\nSaved -> {OUT}/fig_e0_luminosity.[pdf,png]  (+ report/figs/fig_e0_lumi)")


if __name__ == "__main__":
    main()
