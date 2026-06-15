#!/usr/bin/env python3
"""
make_e0_x17_production_fig.py

Third figure in the cross-section series: TOTAL X17 production vs neutron energy,
in the two scenarios set by the X17 quantum numbers.

X17 production = (IPC pair yield) x BR_X17, with BR_X17 = 0.025 (X17 per e+e-
pair, the project benchmark).  Per channel:
    sigma_X17(M1+E1) = sigma_ng * alpha_IPC * BR_X17        (measured sigma_ng)
    sigma_X17(E0)    = sigma_E0 * 1       * BR_X17 = f*sigma_M1*g0+ * BR_X17

Scenarios (the difference is exactly the E0 channel):
  1. VECTOR or SCALAR X17  -> couples to E0 AND M1+E1  (sum both)
  2. PSEUDOSCALAR or AXIAL -> couples to M1+E1 only

Selection rule behind it: a 0+ -> 0+ (E0/monopole) transition can emit a massive
boson only if J^P allows it -- vector(1-)/scalar(0+) yes, pseudoscalar(0-)/
axial(1+) no.  (At the strict per-multipole level the M1/E1 couplings are also
type-dependent; here we take "X17 rides the measured M1+E1 strength" as the
working approximation, so the clean scenario difference is the E0 piece.)

It also computes the CRITICAL numbers: total X17 produced per day in
  - the high-energy window 0.2-2 MeV  (where M1+E1 lives, flash-readout region)
  - the full sub-keV window  E_n < 1 keV  (where E0 lives, thermal-trigger region)
for both scenarios and the f band.

Output: docs/e0_branch/figs/fig_e0_x17_production.{pdf,png}
        (copied to docs/report/figs/fig_e0_x17.{pdf,png})
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
F_E0 = [1e-3, 3e-3, 1e-2]
ALPHA_IPC, ALPHA_TABLE = 3.5e-3, 2.1e-3
BR_X17 = 0.025


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


# ============================ the figure ====================================
def make_figure(En, s_ng):
    boost = g0plus(0.75 * En * 1e-6)
    x17_m1e1 = s_ng * ALPHA_IPC * BR_X17
    x17_e0 = {f: f * s_ng * boost * BR_X17 for f in F_E0}        # x1 pair fraction
    sc1 = {f: x17_m1e1 + x17_e0[f] for f in F_E0}                # vector/scalar: sum

    fig, ax = plt.subplots(figsize=(9.4, 6.0))
    ax.axvspan(1e-3, 1e3, color="C0", alpha=0.07, zorder=0)      # sub-keV window
    ax.axvspan(2e5, 2e6, color="gold", alpha=0.22, zorder=0)     # 0.2-2 MeV window
    ax.text(2e-2, 4e-7, "sub-keV\nwindow", ha="center", fontsize=8.2, color="C0")
    ax.text(6.3e5, 4e-7, "0.2–2 MeV\nwindow", ha="center", fontsize=8.2, color="0.4")

    # scenario 1: total = M1+E1 + E0 (band over f)
    ax.fill_between(En, sc1[F_E0[0]], sc1[F_E0[2]], color="#c0392b", alpha=0.20,
                    label=r"Scenario 1 (vector/scalar): M1+E1 $+$ E0, $f=10^{-3}$–$10^{-2}$")
    ax.loglog(En, sc1[F_E0[1]], color="#c0392b", lw=1.8, ls="--",
              label=r"  …central $f=3\times10^{-3}$")
    # scenario 2: M1+E1 only
    ax.loglog(En, x17_m1e1, color="#2471a3", lw=2.6,
              label=r"Scenario 2 (pseudoscalar/axial): M1+E1 only [measured $\sigma_{n\gamma}$]")

    ax.annotate("E0 adds X17 only here\n(scenario 1, sub-keV)",
                xy=(2e-1, sc1[F_E0[1]][np.argmin(np.abs(En - 2e-1))]),
                xytext=(3e-2, 8e-9), fontsize=8.2, color="#c0392b",
                arrowprops=dict(arrowstyle="->", color="#c0392b"))
    ax.annotate("both scenarios identical here\n(M1+E1, the production peak)",
                xy=(1e6, x17_m1e1[np.argmin(np.abs(En - 1e6))]),
                xytext=(2e3, 3e-9), fontsize=8.2, color="#2471a3",
                arrowprops=dict(arrowstyle="->", color="#2471a3"))

    ax.set_xlabel(r"neutron energy  $E_n$  [eV]")
    ax.set_ylabel(r"X17 production cross section  $\sigma_{\rm X17}$  [barn]")
    ax.set_title("Total X17 production vs neutron energy, by X17 type\n"
                 r"$\sigma_{\rm X17}=$(IPC yield)$\times$BR$_{\rm X17}$;  "
                 "scenarios differ only by the $\\gamma$-dark E0 (monopole) channel")
    ax.set_xlim(1e-3, 1e7)
    ax.set_ylim(1e-10, 1e-6)
    ax.legend(loc="upper center", fontsize=8.0, framealpha=0.95)
    ax.grid(alpha=0.3, which="both")

    OUT.mkdir(parents=True, exist_ok=True)
    REPORT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_e0_x17_production.{ext}", bbox_inches="tight", dpi=150)
        fig.savefig(REPORT / f"fig_e0_x17.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)


# ===================== the window daily statistics ==========================
def window_stats():
    """X17 PRODUCED per day in the 0.2-2 MeV and sub-keV (<1 keV) windows,
    folding the measured per-decade campaign flux (mev_rates.json)."""
    d = json.load(open(RATES))
    ppd = d["pulses_per_day"]
    rows = d["decades"]
    Elo = np.array([r["E_lo_eV"] for r in rows])
    Ehi = np.array([r["E_hi_eV"] for r in rows])
    Ecen = np.sqrt(Elo * Ehi)
    rad_pp = np.array([r["rad_per_pulse"] for r in rows])
    Np = np.array([r["beam_per_pulse"] * r["np_frac_of_beam"] for r in rows])

    _, sng = read_xs("reaction_102"); Eng, _ = read_xs("reaction_102")
    _, snp = read_xs("reaction_103"); Enp, _ = read_xs("reaction_103")
    ratio_th = loginterp(2.53e-8, Eng, sng) / loginterp(2.53e-8, Enp, snp)
    g0 = g0plus(0.75 * Ecen * 1e-6)

    # per-decade X17/day produced
    x17_m1e1_dec = rad_pp * ALPHA_IPC * BR_X17 * ppd
    e0_pairs_dec = Np * F_E0_arr[:, None] * ratio_th * g0 * ppd   # (nf, ndec)
    x17_e0_dec = e0_pairs_dec * BR_X17

    # --- sub-keV (<1 keV): decade-aligned sum ---
    sub = Ehi <= 1e3
    m1e1_sub = x17_m1e1_dec[sub].sum()
    e0_sub = x17_e0_dec[:, sub].sum(axis=1)

    # --- 0.2-2 MeV: M1+E1 from the exact campaign window count ---
    w = d["window_0p2_2MeV"]
    m1e1_hi = w["rad_per_pulse"] * ALPHA_IPC * BR_X17 * ppd
    # E0 in 0.2-2 MeV: log-fraction of the 1e5-1e6 and 1e6-1e7 decades (E0 tiny here)
    i56 = np.argmin(np.abs(Ecen - np.sqrt(1e5 * 1e6)))
    i67 = np.argmin(np.abs(Ecen - np.sqrt(1e6 * 1e7)))
    frac56 = np.log10(1e6 / 2e5) / 1.0      # fraction of 1e5-1e6 in [2e5,1e6]
    frac67 = np.log10(2e6 / 1e6) / 1.0      # fraction of 1e6-1e7 in [1e6,2e6]
    e0_hi = x17_e0_dec[:, i56] * frac56 + x17_e0_dec[:, i67] * frac67

    return dict(ppd=ppd, ratio_th=ratio_th,
                m1e1_sub=m1e1_sub, e0_sub=e0_sub,
                m1e1_hi=m1e1_hi, e0_hi=e0_hi,
                m1e1_hi_table=w["rad_per_pulse"] * ALPHA_TABLE * BR_X17 * ppd)


F_E0_arr = np.array(F_E0)


def main():
    En = np.logspace(-3, 7, 1500)
    _, sng = read_xs("reaction_102"); Eng, _ = read_xs("reaction_102")
    s_ng = loginterp(En, Eng, sng)
    make_figure(En, s_ng)

    st = window_stats()
    print(f"alpha_IPC={ALPHA_IPC}, BR_X17={BR_X17}, pulses/day={st['ppd']:.0f}\n")
    print("=== X17 PRODUCED per day ===")
    print(f"{'window':>16}  {'Scenario 2':>14}  {'Scenario 1 (vector/scalar): M1+E1 + E0(f)':>44}")
    print(f"{'':>16}  {'(M1+E1 only)':>14}  {'f=1e-3':>12} {'f=3e-3':>12} {'f=1e-2':>12}")
    print(f"{'0.2-2 MeV (HE)':>16}  {st['m1e1_hi']:>14.2f}  " +
          "  ".join(f"{st['m1e1_hi']+e:>10.2f}" for e in st['e0_hi']))
    print(f"{'  sub-keV (<1keV)':>16}  {st['m1e1_sub']:>14.4f}  " +
          "  ".join(f"{st['m1e1_sub']+e:>10.4f}" for e in st['e0_sub']))
    print(f"\n  [0.2-2 MeV M1+E1 with table alpha=2.1e-3: {st['m1e1_hi_table']:.1f}/day "
          f"(the committed 22.5/day headline)]")
    print(f"  [E0 alone, sub-keV: f=1e-3..1e-2 -> "
          f"{st['e0_sub'][0]:.3f} .. {st['e0_sub'][2]:.3f} X17/day]")
    print(f"  [E0 alone, 0.2-2 MeV: {st['e0_hi'][0]:.4f} .. {st['e0_hi'][2]:.4f} "
          f"X17/day -- negligible]")
    print("\n  NB: PRODUCED. Recorded ~ x0.196 (MM acceptance); the 0.2-2 MeV window\n"
          "  is in the 10us flash readout, the sub-keV window needs a thermal trigger.")
    print(f"\nSaved -> {OUT}/fig_e0_x17_production.[pdf,png]  (+ report/figs/fig_e0_x17)")


if __name__ == "__main__":
    main()
