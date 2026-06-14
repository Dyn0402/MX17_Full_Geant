#!/usr/bin/env python3
"""
make_e0_pair_yield_fig.py

Folds the 0+ (s-wave, E0) doorway physics into an actual e+e- pair-yield curve
vs neutron energy, and overlays it on the term-1 (M1/E1) curve.  Shows the key
qualitative result: E0 pairs follow the s-wave (n,p)-like distribution and the
sub-threshold 0+ pole -> concentrated at LOW E_n, opposite to term 1 (M1/E1),
which rides the rising sigma_ng/sigma_np and peaks in the MeV region.

Construction (both pieces use measured inputs; only the E0 normalisation f and
the 0+ resonance shape are model):

  term 1 (M1/E1): N1(E) = N_np(E) * [sigma_ng/sigma_np](E) * alpha_IPC
  term 2 (E0):    N0(E) = N_np(E) * f * [sigma_ng/sigma_np]_thermal * g0+(E)

  N_np(E)               : (n,p)t events/pulse/decade from the campaign
                          (= beam_per_pulse * np_frac_of_beam), analysis/mev
  [sigma_ng/sigma_np](E): ENDF/B-VIII.0, data/He3.h5  (rises 1e-8 -> 2e-4)
  alpha_IPC = 3.5e-3    : high-energy IPC coeff (12C/8Be-anchored)
  f = sigma_E0/sigma_M1 : the E0 strength relative to the M1 radiative capture
                          -- THE nuclear unknown; shown for 0.003/0.03/0.3
  g0+(E)                : sub-threshold 0+ Breit-Wigner shape (E_r=-0.368 MeV,
                          Gamma=0.50), normalised to 1 at threshold; makes E0
                          fall toward MeV as the 20.21 pole recedes.

The point is the SHAPE/contrast and the sensitivity to f, not an absolute E0
rate (that needs rho(E0) from 4He(e,e') + the 4He R-matrix).

Output: docs/e0_branch/figs/fig_e0_pair_yield.{pdf,png}
"""
from pathlib import Path
import json
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RATES = ROOT / "analysis" / "mev" / "mev_rates.json"
H5 = ROOT / "data" / "He3.h5"
OUT = ROOT / "docs" / "e0_branch" / "figs"

ALPHA_IPC = 3.5e-3
F_E0 = [0.003, 0.03, 0.3]         # sigma_E0/sigma_M1 illustrative values
F_STYLE = ["-.", "--", "-"]
# sub-threshold 0+ resonance (TUNL): E_x=20.21, threshold 20.578 -> E_r=-0.368
ER_0P, G_0P = 20.21 - 20.578, 0.50
TEMP = "294K"


def endf_ratio(E_MeV):
    """sigma_ng/sigma_np vs energy from ENDF/B-VIII.0."""
    with h5py.File(H5, "r") as f:
        g = f["He3"]
        def xs(rx):
            E = g[f"energy/{TEMP}"][...]
            ds = g[f"reactions/{rx}/{TEMP}/xs"]
            i0 = ds.attrs.get("threshold_idx", 0)
            return E[i0:i0 + ds.shape[0]] * 1e-6, ds[...]
        Eg, ng = xs("reaction_102")
        Ep, npx = xs("reaction_103")
    def il(x, xp, fp):
        return 10 ** np.interp(np.log10(x), np.log10(xp),
                               np.log10(np.clip(fp, 1e-40, None)))
    return il(E_MeV, Eg, ng) / il(E_MeV, Ep, npx)


def g0plus(E_cm):
    num = ER_0P**2 + (G_0P / 2)**2
    return num / ((E_cm - ER_0P)**2 + (G_0P / 2)**2)


def main():
    d = json.load(open(RATES))
    ppd = d["pulses_per_day"]
    rows = d["decades"]
    Elo = np.array([r["E_lo_eV"] for r in rows])
    Ehi = np.array([r["E_hi_eV"] for r in rows])
    Ecen = np.sqrt(Elo * Ehi)                       # eV
    Ecen_MeV = Ecen * 1e-6
    Np = np.array([r["beam_per_pulse"] * r["np_frac_of_beam"] for r in rows])

    ratio = endf_ratio(Ecen_MeV)
    ratio_th = endf_ratio(np.array([2.53e-8]))[0]   # thermal sigma_ng/sigma_np

    term1 = Np * ratio * ALPHA_IPC * ppd            # M1/E1 IPC pairs/day/decade
    Ecm_MeV = 0.75 * Ecen_MeV
    e0 = {f: Np * f * ratio_th * g0plus(Ecm_MeV) * ppd for f in F_E0}

    fig, ax = plt.subplots(figsize=(9, 5.8))
    ax.plot(Ecen, term1, "o-", color="C0", ms=6, lw=2.3,
            label=r"term 1: M1/E1 IPC ($\alpha=3.5\times10^{-3}$, "
                  r"rides $\sigma_{n\gamma}/\sigma_{np}$)")
    for f, ls in zip(F_E0, F_STYLE):
        ax.plot(Ecen, e0[f], ls, color="C3", lw=2,
                label=rf"term 2: E0, $\sigma_{{E0}}/\sigma_{{M1}}={f:g}$")

    ax.axvspan(2e5, 2e6, color="gold", alpha=0.22, zorder=0)
    ax.text(8.5e5, 4e2, "0.2–2 MeV\nwindow", ha="center", fontsize=8.5,
            color="0.4")
    ax.annotate("E0: s-wave + sub-threshold 0$^+$\n→ piles up at LOW $E_n$",
                xy=(3e0, e0[0.03][np.argmin(np.abs(Ecen - 3))]),
                xytext=(2e-2, 2e2), fontsize=9, color="C3",
                arrowprops=dict(arrowstyle="->", color="C3"))
    ax.annotate("M1/E1: rides the rising\n$\\sigma_{n\\gamma}/\\sigma_{np}$ → peaks at MeV",
                xy=(1e6, term1[np.argmin(np.abs(Ecen - 1e6))]),
                xytext=(2e1, 8e-1), fontsize=9, color="C0",
                arrowprops=dict(arrowstyle="->", color="C0"))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"neutron energy  $E_n$  [eV]")
    ax.set_ylabel("e$^+$e$^-$ pairs / day / decade")
    ax.set_title("Where the pairs come from: E0 (low $E_n$) vs M1/E1 (MeV)\n"
                 "E0 normalisation $f=\\sigma_{E0}/\\sigma_{M1}$ is the nuclear "
                 "unknown (needs $\\rho$(E0) + R-matrix)")
    ax.set_ylim(1e-3, 2e3)
    ax.set_xlim(1e-2, 1e8)
    ax.legend(fontsize=8.3, loc="upper left", ncol=1)
    ax.grid(alpha=0.3, which="both")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_e0_pair_yield.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

    # console summary
    print(f"thermal sigma_ng/sigma_np = {ratio_th:.2e}")
    print(f"term 1 (M1/E1) total: {term1.sum():.0f} pairs/day")
    for f in F_E0:
        tot = e0[f].sum()
        # fraction of E0 below 1 keV vs term1 below/above
        lowmask = Ecen < 1e3
        print(f"  E0 (f={f:g}): {tot:6.0f} pairs/day total | "
              f"{e0[f][lowmask].sum():.0f} below 1 keV | "
              f"E0/term1 total = {tot/term1.sum():.2f}")
    print(f"\nSaved -> {OUT}/fig_e0_pair_yield.[pdf,png]")


if __name__ == "__main__":
    main()
