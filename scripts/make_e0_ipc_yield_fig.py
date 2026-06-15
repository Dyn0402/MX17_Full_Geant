#!/usr/bin/env python3
"""
make_e0_ipc_yield_fig.py

Companion to make_e0_cross_section_fig.py: the e+e- pair-production ("IPC
yield") cross section vs neutron energy.  It is simply each capture cross
section times its pair fraction:

    sigma_pair(M1+E1) = sigma_ng(E) x alpha_IPC      (photon channel: only the
                                                      IPC tail makes a pair)
    sigma_pair(E0)    = sigma_E0(E) x 1              (gamma-dark: 100% pairs)

alpha_IPC ~ 3.5e-3 is the high-energy internal-pair coefficient at ~20.6 MeV
(12C/8Be-anchored; band 2.5-4.5e-3; Alberto's table uses 2.1e-3).  sigma_E0 is
f*sigma_M1 with the sub-threshold 0+ boost, f=1e-3..1e-2 (the nuclear unknown).

The point: in the CAPTURE cross section (the other figure) E0 sits 4-5 decades
below M1+E1; multiplying by the pair fraction (x alpha_IPC for M1+E1, x1 for E0)
lifts E0 by 1/alpha_IPC ~ 285 relative to M1+E1, so the PAIR yields are
COMPARABLE at low E_n -- E0 at thermal/sub-keV, M1+E1 at MeV.

Bands: M1+E1 from alpha_IPC (2.5-4.5e-3, ~x1.8); E0 from f (1e-3..1e-2, ~1
decade, the dominant uncertainty).

Output: docs/e0_branch/figs/fig_e0_ipc_yield.{pdf,png}
        (copied to docs/report/figs/fig_e0_ipcyield.{pdf,png})
"""
from pathlib import Path
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
H5 = ROOT / "data" / "He3.h5"
OUT = ROOT / "docs" / "e0_branch" / "figs"
REPORT = ROOT / "docs" / "report" / "figs"
TEMP = "294K"

ER_0P, G_0P = 20.21 - 20.578, 0.50
F_E0 = [1e-3, 3e-3, 1e-2]
ALPHA_IPC, ALPHA_LO, ALPHA_HI = 3.5e-3, 2.5e-3, 4.5e-3


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
    En = np.logspace(-3, 7, 1500)
    Eng, sng = read_xs("reaction_102")
    s_ng = loginterp(En, Eng, sng)
    boost = g0plus(0.75 * En * 1e-6)

    # pair-production cross sections
    pair_m1e1 = s_ng * ALPHA_IPC
    pair_m1e1_lo, pair_m1e1_hi = s_ng * ALPHA_LO, s_ng * ALPHA_HI
    pair_e0 = {f: f * s_ng * boost for f in F_E0}     # x1 (100% pairs)

    fig, ax = plt.subplots(figsize=(9.4, 6.0))
    ax.axvspan(2e5, 2e6, color="gold", alpha=0.20, zorder=0)
    ax.text(6.3e5, 2e-6, "0.2–2 MeV\nwindow", ha="center", fontsize=8.2, color="0.4")

    # M1+E1 pair yield (measured sigma x alpha_IPC), band from alpha
    ax.fill_between(En, pair_m1e1_lo, pair_m1e1_hi, color="#2471a3", alpha=0.20,
                    label=r"M1+E1 band ($\alpha_{\rm IPC}=2.5$–$4.5\times10^{-3}$)")
    ax.loglog(En, pair_m1e1, color="#2471a3", lw=2.6,
              label=r"M1+E1: $\sigma_{n\gamma}\times\alpha_{\rm IPC}$ "
                    r"($\alpha=3.5\times10^{-3}$) [measured $\sigma$]")
    # E0 pair yield (x1), band from f
    ax.fill_between(En, pair_e0[F_E0[0]], pair_e0[F_E0[2]], color="#c0392b",
                    alpha=0.22, label=r"E0 band, $f=10^{-3}$–$10^{-2}$ [ESTIMATE]")
    ax.loglog(En, pair_e0[F_E0[1]], color="#c0392b", lw=1.8, ls="--",
              label=r"E0: $\sigma_{\rm E0}\times1$, $f=3\times10^{-3}$ (central)")

    # thermal comparison marker
    jth = np.argmin(np.abs(En - 0.0253))
    ax.annotate("at thermal the two are COMPARABLE:\n"
                r"E0 (×1) catches the M1+E1 (×$\alpha_{\rm IPC}$)"
                "\nby $1/\\alpha_{\\rm IPC}\\approx285$",
                xy=(En[jth], pair_m1e1[jth]), xytext=(3e-2, 6e-6),
                fontsize=8.0, color="0.25",
                arrowprops=dict(arrowstyle="->", color="0.4"))
    ax.annotate("M1+E1 pairs\npeak at MeV", xy=(2e6, pair_m1e1[np.argmin(np.abs(En-2e6))]),
                xytext=(8e5, 2e-9), fontsize=8.2, color="#2471a3",
                arrowprops=dict(arrowstyle="->", color="#2471a3"))
    ax.annotate("E0 pairs\npile up sub-keV", xy=(3e-1, pair_e0[F_E0[1]][np.argmin(np.abs(En-3e-1))]),
                xytext=(2e-2, 3e-9), fontsize=8.2, color="#c0392b",
                arrowprops=dict(arrowstyle="->", color="#c0392b"))

    ax.set_xlabel(r"neutron energy  $E_n$  [eV]")
    ax.set_ylabel(r"IPC pair-production cross section  $\sigma_{\rm pair}$  [barn]")
    ax.set_title("IPC ($e^+e^-$) yield = capture cross section $\\times$ pair fraction\n"
                 r"M1+E1 $\times\,\alpha_{\rm IPC}$ (measured $\sigma$);  "
                 r"E0 $\times\,1$ ($\gamma$-dark) $\Rightarrow$ comparable yields")
    ax.set_xlim(1e-3, 1e7)
    ax.set_ylim(1e-10, 1e-4)
    ax.legend(loc="upper right", fontsize=8.0, framealpha=0.95)
    ax.grid(alpha=0.3, which="both")

    OUT.mkdir(parents=True, exist_ok=True)
    REPORT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_e0_ipc_yield.{ext}", bbox_inches="tight", dpi=150)
        fig.savefig(REPORT / f"fig_e0_ipcyield.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

    print(f"thermal pair xsec:  M1+E1 = {pair_m1e1[jth]*1e6:.3f} ub "
          f"(band {pair_m1e1_lo[jth]*1e6:.3f}-{pair_m1e1_hi[jth]*1e6:.3f})")
    for f in F_E0:
        print(f"  E0 f={f:<6g}: {pair_e0[f][jth]*1e6:.3f} ub  "
              f"(E0/M1E1 = {pair_e0[f][jth]/pair_m1e1[jth]:.2f})")
    print(f"\nSaved -> {OUT}/fig_e0_ipc_yield.[pdf,png]  (+ report/figs/fig_e0_ipcyield)")


if __name__ == "__main__":
    main()
