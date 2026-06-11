"""Plot 3He(n,p)3H and 3He(n,g)4He cross sections + ratio from ENDF/B-VIII.0.

Default output: ../output/he3_cross_sections.png (standalone use).
Pass an output path to write elsewhere (e.g. the report figs dir):
    python3 scripts/plot_he3_xs.py docs/report/figs/fig_he3_xs.pdf
"""
import os
import sys

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Global parameters ---
_HERE = os.path.dirname(os.path.abspath(__file__))
H5_FILE = os.path.join(_HERE, "../data/He3.h5")  # ENDF/B-VIII.0 (OpenMC/NNDC HDF5)
TEMP = "294K"                  # use 0K grid (room temperature)
E_MIN, E_MAX = 1e-9, 10.0    # MeV
E_THERMAL = 2.53e-8          # MeV
OUT = os.path.join(_HERE, "../output/he3_cross_sections.png")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else OUT
    os.makedirs(os.path.dirname(out), exist_ok=True)
    E, xs_np = load_xs("reaction_103")   # (n,p)
    E2, xs_ng = load_xs("reaction_102")  # (n,gamma)

    # Common grid for ratio
    Eg = np.logspace(np.log10(E_MIN), np.log10(E_MAX), 2000)
    np_i = loglog_interp(Eg, E, xs_np)
    ng_i = loglog_interp(Eg, E2, xs_ng)
    ratio = ng_i / np_i

    s_np_th = loglog_interp(np.array([E_THERMAL]), E, xs_np)[0]
    s_ng_th = loglog_interp(np.array([E_THERMAL]), E2, xs_ng)[0]
    print(f"Thermal (25.3 meV):  sigma(n,p) = {s_np_th:.1f} b")
    print(f"                     sigma(n,g) = {s_ng_th*1e6:.1f} ub")
    print(f"                     ratio (n,g)/(n,p) = {s_ng_th/s_np_th:.2e}")

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8.5, 7.5), sharex=True,
        gridspec_kw={"height_ratios": [1.8, 1], "hspace": 0.07})

    ax1.loglog(Eg, np_i, color="#1f77b4", lw=2,
               label=r"$^3$He(n,p)$^3$H  (MT=103)")
    ax1.loglog(Eg, ng_i, color="#d62728", lw=2,
               label=r"$^3$He(n,$\gamma$)$^4$He  (MT=102)")
    ax1.axvline(E_THERMAL, color="gray", ls=":", lw=1)
    ax1.annotate(f"thermal\n{s_np_th:.0f} b", (E_THERMAL, s_np_th),
                 textcoords="offset points", xytext=(10, 5), fontsize=9)
    ax1.annotate(f"{s_ng_th*1e6:.0f} µb", (E_THERMAL, s_ng_th),
                 textcoords="offset points", xytext=(10, 5), fontsize=9,
                 color="#d62728")
    ax1.set_ylabel("Cross section (barns)")
    ax1.set_title(r"$^3$He + n cross sections (ENDF/B-VIII.0)")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(loc="upper right")

    ax2.loglog(Eg, ratio, color="#2ca02c", lw=2)
    ax2.axvline(E_THERMAL, color="gray", ls=":", lw=1)
    ax2.annotate(f"{s_ng_th/s_np_th:.1e}", (E_THERMAL, s_ng_th / s_np_th),
                 textcoords="offset points", xytext=(10, 5), fontsize=9,
                 color="#2ca02c")
    ax2.set_xlabel("Neutron energy (MeV)")
    ax2.set_ylabel(r"$\sigma_{n\gamma}\,/\,\sigma_{np}$")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.set_xlim(E_MIN, E_MAX)

    # Shade the sub-keV analysis window on both panels
    for ax in (ax1, ax2):
        ax.axvspan(1e-9, 1e-3, color="0.85", alpha=0.45, zorder=0)
    ax2.text(2e-7, 2.5e-8, "sub-keV window", fontsize=8.5, color="0.35")

    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"saved {out}")


def load_xs(rx):
    """Return (E [MeV], sigma [b]) for a reaction group at TEMP."""
    with h5py.File(H5_FILE, "r") as f:
        g = f["He3"]
        E_eV = g[f"energy/{TEMP}"][...]
        ds = g[f"reactions/{rx}/{TEMP}/xs"]
        xs = ds[...]
        i0 = ds.attrs.get("threshold_idx", 0)
    return E_eV[i0:i0 + len(xs)] * 1e-6, xs


def loglog_interp(x, xp, fp):
    fp = np.clip(fp, 1e-30, None)
    return 10 ** np.interp(np.log10(x), np.log10(xp), np.log10(fp))


if __name__ == "__main__":
    main()
