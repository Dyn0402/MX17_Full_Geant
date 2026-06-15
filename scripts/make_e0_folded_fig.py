#!/usr/bin/env python3
"""
make_e0_folded_fig.py

The capstone of the cross-section series: everything folded together into
events per day per decade vs neutron energy, for X17 (left) and IPC e+e- pairs
(right), shown both PRODUCED and RECORDED (x detector acceptance).

Folding (see make_e0_luminosity_fig.py for the bridge):
  M1+E1 :  direct Geant4 (n,gamma) count  rad_per_pulse[decade]
  E0    :  N_(n,p)[decade] * f * (sigma_ng/sigma_np)_th * g0+(E)   (calibrated
           off the breakup; absent from Geant4)
  X17  = (pairs) * BR_X17 ;   IPC pairs = M1+E1*alpha_IPC + E0*1
  recorded = produced * MM-double geometric acceptance
             (0.196 X17, 0.236 IPC; pairs_v2 at-rest study).

Two scenarios for X17 by quantum numbers:
  scenario 2 (pseudoscalar/axial): M1+E1 only
  scenario 1 (vector/scalar):      M1+E1 + E0   (band over f=1e-3..1e-2)

Text boxes give the per-day SUMS in the 0.2-2 MeV (high-energy) and sub-keV
(<1 keV) windows, produced and recorded.

Output: docs/e0_branch/figs/fig_e0_folded.{pdf,png}
        (copied to docs/report/figs/fig_e0_folded.{pdf,png})
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
F_E0 = [1e-3, 3e-3, 1e-2]
ACC_X17, ACC_IPC = 0.196, 0.236


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
    ppd = d["pulses_per_day"]
    rows = d["decades"]
    Elo = np.array([r["E_lo_eV"] for r in rows])
    Ehi = np.array([r["E_hi_eV"] for r in rows])
    Ecen = np.sqrt(Elo * Ehi)
    rad_pp = np.array([r["rad_per_pulse"] for r in rows])
    Np = np.array([r["beam_per_pulse"] * r["np_frac_of_beam"] for r in rows])

    Eng, sng = read_xs("reaction_102")
    Enp, snp = read_xs("reaction_103")
    ratio_th = loginterp(2.53e-8, Eng, sng) / loginterp(2.53e-8, Enp, snp)
    g0 = g0plus(0.75 * Ecen * 1e-6)

    # per-decade IPC pairs/day
    ipc_m1e1 = rad_pp * ALPHA_IPC * ppd
    ipc_e0 = {f: Np * f * ratio_th * g0 * ppd for f in F_E0}
    # X17 = pairs * BR
    x17_m1e1 = ipc_m1e1 * BR_X17
    x17_e0 = {f: ipc_e0[f] * BR_X17 for f in F_E0}

    # ---- window sums -------------------------------------------------------
    sub = Ehi <= 1e3
    w = d["window_0p2_2MeV"]
    frac56, frac67 = np.log10(1e6 / 2e5), np.log10(2e6 / 1e6)
    i56 = np.argmin(np.abs(Ecen - np.sqrt(1e5 * 1e6)))
    i67 = np.argmin(np.abs(Ecen - np.sqrt(1e6 * 1e7)))

    def sums(m1e1_dec, e0_dec, acc):
        # M1+E1 high-E: exact campaign window count x the same per-(n,g) factor
        # (= alpha_IPC*[BR or 1]) that built the per-decade curve.
        per_rad = m1e1_dec[i67] / (rad_pp[i67] * ppd)        # = alpha * (BR or 1)
        m1e1_hi = w["rad_per_pulse"] * per_rad * ppd
        m1e1_sub = m1e1_dec[sub].sum()
        e0_hi = e0_dec[:, i56] * frac56 + e0_dec[:, i67] * frac67
        e0_sub = e0_dec[:, sub].sum(axis=1)
        return dict(m1e1_hi=m1e1_hi, m1e1_sub=m1e1_sub, e0_hi=e0_hi, e0_sub=e0_sub,
                    acc=acc)

    e0_x17_arr = np.array([x17_e0[f] for f in F_E0])
    e0_ipc_arr = np.array([ipc_e0[f] for f in F_E0])
    sx = sums(x17_m1e1, e0_x17_arr, ACC_X17)
    si = sums(ipc_m1e1, e0_ipc_arr, ACC_IPC)

    # ======================= figure ========================================
    fig, (axX, axI) = plt.subplots(1, 2, figsize=(13.6, 6.2), sharex=True, sharey=False)

    def panel(ax, m1e1, e0_arr, acc, title, ylab):
        ax.axvspan(1e-3, 1e3, color="C0", alpha=0.06, zorder=0)
        ax.axvspan(2e5, 2e6, color="gold", alpha=0.20, zorder=0)
        sc1_lo, sc1_hi = m1e1 + e0_arr[0], m1e1 + e0_arr[2]
        sc1_mid = m1e1 + e0_arr[1]
        # produced
        ax.fill_between(Ecen, sc1_lo, sc1_hi, color="#c0392b", alpha=0.20,
                        label="scenario 1 (vector/scalar): +E0, $f=10^{-3}$–$10^{-2}$")
        ax.plot(Ecen, sc1_mid, "--", color="#c0392b", lw=1.6)
        ax.plot(Ecen, m1e1, "o-", color="#2471a3", ms=5, lw=2.2,
                label="scenario 2 (pseudo/axial): M1+E1 only [produced]")
        # recorded (x acceptance)
        ax.plot(Ecen, m1e1 * acc, "s--", color="#2471a3", ms=4, lw=1.4, alpha=0.8,
                label=f"recorded = ×{acc} (MM acceptance)")
        ax.fill_between(Ecen, sc1_lo * acc, sc1_hi * acc, color="#c0392b",
                        alpha=0.12, hatch="//")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"neutron energy  $E_n$  [eV]")
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.set_xlim(1e-3, 1e7)
        ax.grid(alpha=0.3, which="both")
        ax.legend(loc="lower center", fontsize=7.6, framealpha=0.95)

    panel(axX, x17_m1e1, e0_x17_arr, ACC_X17,
          "X17 production", "X17 / day / decade")
    panel(axI, ipc_m1e1, e0_ipc_arr, ACC_IPC,
          "IPC ($e^+e^-$) pairs", "IPC pairs / day / decade")

    # ---- text boxes with the window sums ----------------------------------
    def box(ax, S, label, unit):
        hi_p = S['m1e1_hi'];  hi_r = hi_p * S['acc']
        sub2_p = S['m1e1_sub']; sub2_r = sub2_p * S['acc']
        sub1_lo_p = S['m1e1_sub'] + S['e0_sub'][0]
        sub1_hi_p = S['m1e1_sub'] + S['e0_sub'][2]
        txt = (f"{label} / day   (produced → recorded ×{S['acc']})\n"
               f"  0.2–2 MeV : {hi_p:,.0f} → {hi_r:,.0f}   (both scen.)\n"
               f"  sub-keV, scen.2 : {sub2_p:.2f} → {sub2_r:.2f}\n"
               f"  sub-keV, scen.1 : {sub1_lo_p:.2f}–{sub1_hi_p:.2f}"
               f" → {sub1_lo_p*S['acc']:.2f}–{sub1_hi_p*S['acc']:.2f}")
        ax.text(0.025, 0.975, txt, transform=ax.transAxes, va="top", ha="left",
                fontsize=7.8, family="monospace",
                bbox=dict(boxstyle="round", fc="#fffbe6", ec="0.5", alpha=0.95))

    box(axX, sx, "X17", "X17")
    box(axI, si, "IPC pairs", "pairs")

    axX.set_ylim(1e-3, 1e2)
    axI.set_ylim(1e-1, 1e4)

    fig.suptitle("Folded result: X17 and IPC per day per decade, produced and "
                 "recorded (×acceptance)\n"
                 "gold = 0.2–2 MeV window, blue = sub-keV; sub-keV recording also "
                 "needs a thermal trigger (not applied here)", fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    OUT.mkdir(parents=True, exist_ok=True)
    REPORT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_e0_folded.{ext}", bbox_inches="tight", dpi=150)
        fig.savefig(REPORT / f"fig_e0_folded.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

    # console
    for name, S in [("X17", sx), ("IPC", si)]:
        print(f"=== {name} (acc={S['acc']}) ===")
        print(f"  0.2-2 MeV : produced {S['m1e1_hi']:.2f}  recorded {S['m1e1_hi']*S['acc']:.2f}")
        print(f"  sub-keV scen.2: produced {S['m1e1_sub']:.3f}  recorded {S['m1e1_sub']*S['acc']:.3f}")
        print(f"  sub-keV scen.1: produced {S['m1e1_sub']+S['e0_sub'][0]:.3f}"
              f"–{S['m1e1_sub']+S['e0_sub'][2]:.3f}  recorded "
              f"{(S['m1e1_sub']+S['e0_sub'][0])*S['acc']:.3f}"
              f"–{(S['m1e1_sub']+S['e0_sub'][2])*S['acc']:.3f}")
    print(f"\nSaved -> {OUT}/fig_e0_folded.[pdf,png]  (+ report/figs/fig_e0_folded)")


if __name__ == "__main__":
    main()
