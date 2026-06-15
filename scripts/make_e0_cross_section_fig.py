#!/usr/bin/env python3
"""
make_e0_cross_section_fig.py

The "final cross section" plot for the E0 thread: the radiative-capture cross
sections of n+3He vs neutron energy, separating what is MEASURED from what is
ESTIMATED, with the E0 channel shown as an uncertainty band.

Physics (see doorway_states_note.md): the capture cross section factorises as
    sigma(n,gamma) = sigma_form(Jpi) x [Gamma_gamma / Gamma_tot]   (form x radiate)
and per channel sigma_X = sigma_form x Gamma_X/Gamma_tot.  We plot:

  - sigma_np  : 3He(n,p)t -- the dominant (breakup) channel, ENDF (context).
  - sigma_ng  : 3He(n,gamma)4He total RADIATIVE capture (M1+E1), ENDF/B-VIII.0
                = the photon-emitting capture.  MEASURED/evaluated (thermal
                = 55 ub, the classic value).  Shows the M1-dominated 1/v at low
                E and the E1/GDR rise toward MeV.  This is the KNOWN curve.
  - sigma_E0  : the gamma-dark E0 capture = f * sigma_M1, with
                f = sigma_E0/sigma_M1 the nuclear unknown (band 1e-3..1e-2) and
                the sub-threshold 0+ boost g0+(E) steepening it toward low E.
                sigma_M1 ~ sigma_ng at low E (M1-dominated), where E0 lives.
                This is the ESTIMATED curve -- shown as a band.

The band is the f uncertainty (the dominant one, ~1 decade); the shape (1/v x
sub-threshold boost) is robust.  f itself needs a 4He R-matrix + (e,e') monopole
normalisation to pin (roadmap Phase 2).

Output: docs/e0_branch/figs/fig_e0_cross_section.{pdf,png}
        (copied to docs/report/figs/fig_e0_xsec.{pdf,png})
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

# sub-threshold 0+ (TUNL): E_x=20.21, n+3He threshold 20.578 -> E_r=-0.368 MeV
ER_0P, G_0P = 20.21 - 20.578, 0.50
F_E0 = [1e-3, 3e-3, 1e-2]      # sigma_E0/sigma_M1 bracket (lo, mid, hi)


def read_xs(rx):
    with h5py.File(H5, "r") as f:
        g = f["He3"]
        E = g[f"energy/{TEMP}"][...]
        ds = g[f"reactions/{rx}/{TEMP}/xs"]
        i0 = ds.attrs.get("threshold_idx", 0)
        xs = ds[...]
        return E[i0:i0 + xs.shape[0]], xs        # eV, barns


def loginterp(x, xp, fp):
    return 10 ** np.interp(np.log10(x), np.log10(xp),
                           np.log10(np.clip(fp, 1e-40, None)))


def g0plus(E_cm_MeV):
    """sub-threshold 0+ Breit-Wigner boost, normalised to 1 at E_cm=0 (thermal)."""
    num = ER_0P**2 + (G_0P / 2)**2
    return num / ((E_cm_MeV - ER_0P)**2 + (G_0P / 2)**2)


def main():
    En = np.logspace(-3, 7, 1500)           # eV, thermal -> 10 MeV
    Eng, sng = read_xs("reaction_102")      # (n,gamma)
    Enp, snp = read_xs("reaction_103")      # (n,p)
    s_ng = loginterp(En, Eng, sng)
    s_np = loginterp(En, Enp, snp)

    # E0: sigma_E0 = f * sigma_M1 * g0+(E); sigma_M1 ~ sigma_ng at low E.
    Ecm_MeV = 0.75 * En * 1e-6
    boost = g0plus(Ecm_MeV)
    s_e0 = {f: f * s_ng * boost for f in F_E0}

    fig, ax = plt.subplots(figsize=(9.4, 6.0))
    ax.axvspan(2e5, 2e6, color="gold", alpha=0.20, zorder=0)
    ax.text(6.3e5, 3e3, "0.2–2 MeV\nwindow", ha="center", fontsize=8.2, color="0.4")

    ax.loglog(En, s_np, color="0.6", lw=1.8, ls=":",
              label=r"$\sigma_{np}$ (n,p)t — dominant breakup [ENDF]")
    ax.loglog(En, s_ng, color="#2471a3", lw=2.6,
              label=r"$\sigma_{n\gamma}$ radiative capture, M1+E1 [ENDF, MEASURED]")
    ax.fill_between(En, s_e0[F_E0[0]], s_e0[F_E0[2]], color="#c0392b", alpha=0.22,
                    label=r"$\sigma_{\rm E0}=f\,\sigma_{\rm M1}$, $f=10^{-3}$–$10^{-2}$ "
                          r"[ESTIMATE, $\gamma$-dark]")
    ax.loglog(En, s_e0[F_E0[1]], color="#c0392b", lw=1.8, ls="--",
              label=r"$\sigma_{\rm E0}$, $f=3\times10^{-3}$ (central)")

    # annotations: the M1 -> E1 story on the measured curve
    ax.annotate("M1-dominated,\n$1/v$ (direct)", xy=(3e0, loginterp(3e0, Eng, sng)),
                xytext=(3e-3, 3e-3), fontsize=8.2, color="#2471a3",
                arrowprops=dict(arrowstyle="->", color="#2471a3"))
    ax.annotate("E1 / GDR\nturning on", xy=(3e6, loginterp(3e6, Eng, sng)),
                xytext=(8e5, 8e-6), fontsize=8.2, color="#2471a3",
                arrowprops=dict(arrowstyle="->", color="#2471a3"))
    ax.annotate("sub-threshold $0^+$ boost\n$\\Rightarrow$ E0 steepens at low $E_n$",
                xy=(1e0, s_e0[F_E0[1]][np.argmin(np.abs(En - 1e0))]),
                xytext=(2e-2, 2e-9), fontsize=8.2, color="#c0392b",
                arrowprops=dict(arrowstyle="->", color="#c0392b"))

    # thermal markers
    jth = np.argmin(np.abs(En - 0.0253))
    ax.plot(En[jth], s_ng[jth], "o", color="#2471a3", ms=5)
    ax.text(En[jth] * 1.5, s_ng[jth] * 0.32, f"{s_ng[jth]*1e6:.0f} µb (thermal)",
            fontsize=7.6, color="#2471a3", ha="left", va="top")

    ax.set_xlabel(r"neutron energy  $E_n$  [eV]")
    ax.set_ylabel(r"cross section  $\sigma$  [barn]")
    ax.set_title("n+$^3$He cross sections: what is measured ($\\sigma_{n\\gamma}$) "
                 "vs what is estimated ($\\sigma_{\\rm E0}$)\n"
                 "$\\sigma_X=\\sigma_{\\rm form}\\times\\Gamma_X/\\Gamma_{\\rm tot}$;"
                 " the E0 band is the $f$ unknown ($\\sim$1 decade)")
    ax.set_xlim(1e-3, 1e7)
    ax.set_ylim(1e-10, 1e5)
    ax.legend(loc="upper right", fontsize=8.0, framealpha=0.95)
    ax.grid(alpha=0.3, which="both")

    OUT.mkdir(parents=True, exist_ok=True)
    REPORT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_e0_cross_section.{ext}", bbox_inches="tight", dpi=150)
        fig.savefig(REPORT / f"fig_e0_xsec.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

    # console
    print(f"thermal:  sigma_np = {s_np[jth]:.4g} b   sigma_ng = {s_ng[jth]*1e6:.1f} ub"
          f"   ratio ng/np = {s_ng[jth]/s_np[jth]:.2e}")
    for f in F_E0:
        print(f"  f={f:<6g}: sigma_E0(thermal) = {s_e0[f][jth]*1e6:.3g} ub "
              f"= {s_e0[f][jth]/s_np[jth]:.2e} x sigma_np")
    print(f"\nSaved -> {OUT}/fig_e0_cross_section.[pdf,png]  (+ report/figs/fig_e0_xsec)")


if __name__ == "__main__":
    main()
