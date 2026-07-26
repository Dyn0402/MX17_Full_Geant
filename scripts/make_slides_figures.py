#!/usr/bin/env python3
"""make_slides_figures.py — figures for the INFN/collaboration slide deck.

Reuses the validated pair kinematics + measured-response smearing machinery
from scripts/make_angular_resolution_figs.py, and the directly-counted MeV
capture energies from analysis/mev/mev_captures.npz + mev_rates.json.

Two run scenarios (decided 2026-06-16):
  * JULY  (current hardware, optimistic): neutron window 0.2-0.7 MeV, set by a
           CONSERVATIVE gamma-flash recovery of ~1.7 us at L = 19.5 m.
  * LS3   (post hardware optimisation):   window extended to 0.2-2.0 MeV
           (recovery ~1.0 us) -> reach to 2 MeV neutrons, ~3.4x the statistics.

Normalisation: alpha_IPC = 3.5e-3 (12C/8Be-anchored), 30-day run, MM-double
acceptance NOSE-FIRST 27.8% (X17) / 27.4% (IPC) (2026-07-26, analysis/pairs_nose;
was 19.6/23.6 in the 2026-06-18 pre-final-geometry projection).  Recorded yields
computed here from the raw 714 He3(n,g) energies, so any window is exact.
Asimov Z (stat-only): July 3.5 sigma, LS3 6.4 sigma (was 2.6/4.9); the higher
acceptance lifts significance x1.31.  IPC-shape systematic still sets the CL.

Outputs -> docs/slides/figs/ :
  fig_stacked_july.png      stacked smeared e+e- opening-angle spectrum, July
  fig_stacked_ls3.png       same, post-LS3 window
  fig_stacked_compare.png   July | LS3 side by side (deck centrepiece)
  fig_production_windows.png X17/day vs neutron energy, both windows shaded

Usage:  /home/dylan/PycharmProjects/nTof_x17/venv/bin/python \
            scripts/make_slides_figures.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from make_angular_resolution_figs import (          # noqa: E402
    M_X17, sample_pairs, sample_ipc_masses,
    make_psi_sampler, smear_pair_response,
)

OUT = REPO / "docs/slides/figs"
NPZ = REPO / "analysis/mev/mev_captures.npz"
RATES = REPO / "analysis/mev/mev_rates.json"
RESP = REPO / "analysis/pairs_v2/geant4_response.json"

# ── physics / run constants ──────────────────────────────────────────────────
ALPHA_IPC = 3.5e-3        # 12C/8Be-anchored IPC coefficient (default)
BR_X17 = 0.025            # X17 / IPC branching
ACC_X17, ACC_IPC = 0.278, 0.274   # MM-double acceptance, NOSE-FIRST final geom
                                  # (analysis/pairs_nose; was 0.196/0.236 pre-final)
DAYS = 30
M_N_EV = 939.565e6        # neutron mass [eV/c^2]
L_M = 19.5                # EAR2 flight path [m]
C_MS = 2.99792458e8

C_IPC = "#9bb8d4"
C_X17 = "#e84040"
C_TRUTH = "#2ca02c"


def tof_us(E_eV):
    """Non-relativistic neutron time of flight to L = 19.5 m, in us."""
    return L_M / (C_MS * np.sqrt(2.0 * np.asarray(E_eV) / M_N_EV)) * 1e6


def window_yields(E, w, lo, hi):
    """Recorded 30-day X17 / IPC yields for a neutron window [lo, hi] eV."""
    n_rad = int(((E >= lo) & (E < hi)).sum())
    rad_pp = n_rad * w                                   # captures / pulse
    ipc_pp = rad_pp * ALPHA_IPC
    x17_pp = ipc_pp * BR_X17
    ppd = JSON_RATES["pulses_per_day"]
    return dict(
        n_rad=n_rad, rad_pp=rad_pp,
        x17_prod_day=x17_pp * ppd, ipc_prod_day=ipc_pp * ppd,
        x17_rec=x17_pp * ppd * ACC_X17 * DAYS,
        ipc_rec=ipc_pp * ppd * ACC_IPC * DAYS,
        tof_lo=float(tof_us(hi)), tof_hi=float(tof_us(lo)),
    )


# ── the stacked smeared spectrum (one scenario) ──────────────────────────────
#   Base binning is 90 uniform 2° bins over 0-180°.  `rebin` groups that many
#   adjacent base bins into one (rebin=4 -> 8° bins): 90 is not divisible by 4,
#   so the final coarse bin (176-180°) holds the 2 leftover base bins; the
#   opening-angle distribution is empty there, so this edge effect is cosmetic.
BASE_BIN_DEG = 2.0


def rebinned_edges(rebin):
    """Bin edges for the requested rebin factor (1 -> the native 2° grid)."""
    fine = np.linspace(0, 180, 91)
    if rebin <= 1:
        return fine
    return np.unique(np.append(fine[::rebin], 180.0))


def stacked_panel(ax, d, n_x17, n_ipc, n=2_000_000, seed=44, show_truth=True,
                  rebin=1):
    """Draw IPC + X17 stacked, both smeared by the measured P(psi|KE) response
    (best estimator = target-centre chord), normalised to recorded yields.
    `rebin` coarsens the 2° base binning by that integer factor."""
    rng = np.random.default_rng(seed)
    sampler = make_psi_sampler(d, "nomline")
    th_x, ke1x, ke2x, u1x, u2x = sample_pairs(M_X17, n, rng)
    mi = sample_ipc_masses(n, rng)
    th_i, ke1i, ke2i, u1i, u2i = sample_pairs(mi, n, rng)
    th_x_s = smear_pair_response(ke1x, ke2x, u1x, u2x, sampler, rng)
    th_i_s = smear_pair_response(ke1i, ke2i, u1i, u2i, sampler, rng)

    bins = rebinned_edges(rebin)
    cen = 0.5 * (bins[:-1] + bins[1:])

    def counts(th, total):
        h, _ = np.histogram(th, bins=bins)
        return h / h.sum() * total

    b = counts(th_i_s, n_ipc)
    s = counts(th_x_s, n_x17)
    ax.fill_between(cen, 0, b, step="mid", color=C_IPC, alpha=0.95,
                    label=f"IPC background  ({n_ipc:.0f})")
    ax.fill_between(cen, b, b + s, step="mid", color=C_X17, alpha=0.95,
                    label=f"X17 signal  ({n_x17:.0f})")
    ax.step(cen, b + s, where="mid", color="k", lw=1.0)
    if show_truth:
        # undiluted X17 shoulder (truth), riding on the IPC baseline, to show
        # what the capsule multiple-scattering resolution costs the peak.
        st = counts(th_x, n_x17)
        ax.step(cen, b + st, where="mid", color=C_TRUTH, lw=1.6, ls=":",
                label="X17 truth shoulder (no smearing)")
    ax.set_xlim(0, 180)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("e$^+$e$^-$ opening angle  [deg]")
    ax.grid(True, alpha=0.3)
    return cen, b, s


def fig_stacked(scenario, d, E, w, lo, hi, title_extra="", rebin=1):
    y = window_yields(E, w, lo, hi)
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    stacked_panel(ax, d, y["x17_rec"], y["ipc_rec"], rebin=rebin)
    bw = BASE_BIN_DEG * rebin
    ax.set_ylabel(f"recorded events / {bw:g}$^\\circ$  (30-day run)")
    ax.set_title(
        f"{scenario}: stacked e$^+$e$^-$ opening-angle spectrum{title_extra}\n"
        f"$E_n$ = {lo/1e6:g}-{hi/1e6:g} MeV  (TOF {y['tof_lo']:.1f}-{y['tof_hi']:.1f} "
        f"$\\mu$s)   |   {y['x17_rec']:.0f} X17 on {y['ipc_rec']:.0f} IPC, "
        f"S/B = {y['x17_rec']/y['ipc_rec']:.3f}",
        fontsize=10.5)
    ax.legend(fontsize=10, loc="upper right")
    fig.tight_layout()
    suffix = f"_rebin{rebin}" if rebin > 1 else ""
    fig.savefig(OUT / f"fig_stacked_{scenario.lower().split()[0]}{suffix}.png",
                dpi=160)
    plt.close(fig)
    return y


def fig_compare(d, E, w, rebin=1):
    yj = window_yields(E, w, 2e5, 7e5)
    yl = window_yields(E, w, 2e5, 2e6)
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.7), sharey=False)
    for ax, y, lab, lo, hi in [
        (axes[0], yj, "July (current hardware, optimistic)", 2e5, 7e5),
        (axes[1], yl, "Post-LS3 (reach extended to 2 MeV)", 2e5, 2e6),
    ]:
        stacked_panel(ax, d, y["x17_rec"], y["ipc_rec"], rebin=rebin)
        ax.set_title(
            f"{lab}\n$E_n$ {lo/1e6:g}-{hi/1e6:g} MeV "
            f"(TOF {y['tof_lo']:.1f}-{y['tof_hi']:.1f} $\\mu$s)  |  "
            f"{y['x17_rec']:.0f} X17 / {y['ipc_rec']:.0f} IPC", fontsize=10.5)
        ax.legend(fontsize=9.5, loc="upper right")
    bw = BASE_BIN_DEG * rebin
    axes[0].set_ylabel(f"recorded events / {bw:g}$^\\circ$  (30-day run)")
    fig.suptitle(
        "Expected recorded opening-angle spectrum: July vs post-LS3   "
        "($\\alpha_{IPC}=3.5\\times10^{-3}$, MM-double acceptance, "
        "best-estimator smearing)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    suffix = f"_rebin{rebin}" if rebin > 1 else ""
    fig.savefig(OUT / f"fig_stacked_compare{suffix}.png", dpi=160)
    plt.close(fig)
    return yj, yl


# ── production vs neutron energy, with both windows shaded ────────────────────
def fig_production_windows(E, w):
    ppd = JSON_RATES["pulses_per_day"]
    k = ALPHA_IPC * BR_X17 * ppd                 # captures/pulse -> X17/day
    lo_e, hi_e = 1e3, 5e7                         # plotted range [eV]
    edges = np.logspace(np.log10(lo_e), np.log10(hi_e), 26)   # 5 bins/decade
    counts, _ = np.histogram(E, bins=edges)
    x17_day = counts * w * k
    widths = np.diff(edges)
    rec_july = ((E >= 2e5) & (E < 7e5)).sum() * w * k * ACC_X17 * DAYS
    rec_ls3 = ((E >= 2e5) & (E < 2e6)).sum() * w * k * ACC_X17 * DAYS

    fig, ax = plt.subplots(figsize=(9.4, 5.4))
    ax.bar(edges[:-1], x17_day, width=widths, align="edge", color="#9aa0a6",
           edgecolor="k", lw=0.4, label="X17 produced / day (per bin)")
    ax.axvspan(2e5, 7e5, color="#2ca02c", alpha=0.20, zorder=0,
               label="July window  0.2-0.7 MeV")
    ax.axvspan(7e5, 2e6, color="#1f77b4", alpha=0.20, zorder=0,
               label="LS3 extension  0.7-2 MeV")
    for x, txt in [(7e5, "July edge\n~0.7 MeV"), (2e6, "LS3 edge\n2 MeV")]:
        ax.axvline(x, color="k", lw=1.0, ls="--")
    ax.set_xscale("log")
    ax.set_xlim(lo_e, hi_e)
    ymax = x17_day.max() * 1.25
    ax.set_ylim(0, ymax)
    ax.set_xlabel("neutron energy  $E_n$  [eV]")
    ax.set_ylabel("X17 produced / day  ($\\alpha_{IPC}=3.5\\times10^{-3}$)")
    ax.set_title("Where the X17 rate lives: direct $^3$He(n,$\\gamma$) counts "
                 "vs neutron energy\nLS3 extends the $\\gamma$-flash-limited "
                 "reach ~0.7 MeV $\\rightarrow$ 2 MeV", fontsize=10.5)
    ax.text(0.02, 0.7,
            f"recorded X17 / 30 d\n  July  0.2-0.7 MeV : {rec_july:.0f}\n"
            f"  LS3   0.2-2 MeV : {rec_ls3:.0f}   ($\\times${rec_ls3/rec_july:.1f})",
            transform=ax.transAxes, ha="left", va="top", fontsize=10,
            family="monospace",
            bbox=dict(boxstyle="round", fc="white", ec="#888"))
    secax = ax.secondary_xaxis(
        "top", functions=(lambda e: tof_us(np.clip(e, 1.0, None)),
                          lambda t: M_N_EV / 2 *
                          (L_M / (C_MS * np.clip(t, 1e-6, None) * 1e-6)) ** 2))
    secax.set_xlabel("neutron TOF to EAR2 (L = 19.5 m)  [$\\mu$s]")
    ax.legend(fontsize=9.5, loc="upper left")
    ax.grid(True, alpha=0.3, which="major")
    fig.tight_layout()
    fig.savefig(OUT / "fig_production_windows.png", dpi=160)
    plt.close(fig)


def main():
    global JSON_RATES
    OUT.mkdir(parents=True, exist_ok=True)
    JSON_RATES = json.load(open(RATES))
    d = json.load(open(RESP))
    npz = np.load(NPZ)
    E = npz["rad_events"].ravel()
    w = JSON_RATES["weight_per_neutron"]

    yj = fig_stacked("July", d, E, w, 2e5, 7e5)
    yl = fig_stacked("LS3", d, E, w, 2e5, 2e6)
    fig_compare(d, E, w)
    fig_production_windows(E, w)

    # rebin-by-4 versions (8° bins) for the backup slides (colleague request)
    fig_stacked("July", d, E, w, 2e5, 7e5, rebin=4)
    fig_stacked("LS3", d, E, w, 2e5, 2e6, rebin=4)
    fig_compare(d, E, w, rebin=4)

    print("Recorded yields (30-day run, alpha_IPC = 3.5e-3):")
    for nm, y in [("July 0.2-0.7 MeV", yj), ("LS3 0.2-2 MeV", yl)]:
        print(f"  {nm:18s}  X17={y['x17_rec']:6.1f}  IPC={y['ipc_rec']:7.0f}  "
              f"S/B={y['x17_rec']/y['ipc_rec']:.3f}  "
              f"(N_capt={y['n_rad']}, TOF {y['tof_lo']:.2f}-{y['tof_hi']:.2f} us)")

    # Asimov significance (stat-only; IPC-shape systematic NOT included and
    # dominates the achievable CL) from the binned smeared opening-angle templates
    def asimov_Z(s, b):
        m = b > 0
        return float(np.sqrt(2.0 * np.sum((s[m] + b[m]) * np.log1p(s[m] / b[m]) - s[m])))
    _fig, _ax = plt.subplots()
    print("Asimov Z (stat-only, binned profile-likelihood over the smeared templates):")
    for nm, lo, hi in [("July 0.2-0.7 MeV", 2e5, 7e5), ("LS3 0.2-2 MeV", 2e5, 2e6)]:
        y = window_yields(E, w, lo, hi)
        _ax.clear()
        _, b_, s_ = stacked_panel(_ax, d, y["x17_rec"], y["ipc_rec"])
        print(f"  {nm:18s}  Z = {asimov_Z(s_, b_):.2f} sigma")
    plt.close(_fig)
    print(f"figures written to {OUT}")


if __name__ == "__main__":
    main()
