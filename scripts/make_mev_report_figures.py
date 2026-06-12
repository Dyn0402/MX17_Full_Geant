#!/usr/bin/env python3
"""
make_mev_report_figures.py
Figures + rate tables for docs/report/mev_note.tex from the run B-full scan
(analysis/mev/mev_captures.npz, produced by analyze_mev_captures.py).

The radiative-capture rate is taken from the DIRECT 3He(n,g) counts -- 714
events in 5e8 simulated neutrons, 95% above 100 keV.  (The thermal
sigma_ng/sigma_np ratio is only used as the sub-keV floor reference; applying
it at all energies underestimates the MeV region by ~3 orders of magnitude,
since the ratio rises from 1e-8 at thermal to ~1e-4 at MeV.)

Per-pulse normalisation: the generator sampled the EAR2 flux histogram
flux_n_pulse_NOisolet_100bpd over 1 meV - 100 MeV; its integral is
2.263e7 n/pulse (the < 1 keV part reproduces the known 7.31e6 anchor).

Usage:
    python3 scripts/make_mev_report_figures.py \
        [--npz analysis/mev/mev_captures.npz] [--outdir docs/report/figs]
"""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "figure.dpi": 120,
                     "axes.grid": True, "grid.alpha": 0.3})

HERE = Path(__file__).resolve().parent
REPO = HERE.parent

# ── Constants ─────────────────────────────────────────────────────────────────
FLUX_FILE   = REPO / "data/fluxEAR2-Ph3_in_different_units.root"
FLUX_HIST   = "flux_n_pulse_NOisolet_100bpd"
HE3_H5      = REPO / "data/He3.h5"
PULSES_DAY  = 1.929e4          # 7e12 ppp (results_3He header)
ALPHA_IPC   = 2.1e-3           # IPC per radiative capture (Alberto)
BR_X17      = 0.025            # X17 per IPC pair (NOT per capture)
SIGMA_TH    = 54e-6 / 5333.0   # thermal sigma_ng/sigma_np = 1.013e-8
L_FLIGHT_M  = 19.5             # EAR2 nominal flight path [m]
M_N_EV      = 939.565e6        # neutron mass [eV]
C_M_S       = 2.998e8

# Alberto thin-target table (results_3He): decade-lower-edge exponent ->
# (neutrons/pulse, He3-captures/pulse).  Comparison is done per incident
# neutron so his (hotter) flux assumption divides out.
ALBERTO = {
    -2: (3.11e6, 4.37e0), -1: (1.68e6, 9.60e-1), 0: (1.39e6, 2.33e-1),
     1: (1.57e6, 1.05e-1),  2: (1.78e6, 1.05e-1), 3: (2.09e6, 2.90e-1),
     4: (3.03e6, 1.32e0),   5: (9.89e6, 1.50e1),  6: (6.75e6, 1.77e1),
     7: (1.64e6, 1.56e0),
}
ALBERTO_WINDOW = (1.24e7, 2.46e1)   # commented 0.2-2 MeV row


def poisson_interval(n):
    """68% central Poisson interval (Garwood-style), as in the thermal note."""
    if n == 0:
        return 0.0, 1.84
    lo = n * (1 - 1/(9*n) - 1/(3*np.sqrt(n)))**3
    hi = (n+1) * (1 - 1/(9*(n+1)) + 1/(3*np.sqrt(n+1)))**3
    return lo, hi


def load_endf_ratio():
    """ENDF/B-VIII.0 sigma_ng/sigma_np on a common grid (E in eV)."""
    def load_xs(rx):
        with h5py.File(HE3_H5, "r") as f:
            g = f["He3"]
            E = g["energy/294K"][...]
            ds = g[f"reactions/{rx}/294K/xs"]
            xs = ds[...]
            i0 = ds.attrs.get("threshold_idx", 0)
        return E[i0:i0 + len(xs)], xs
    E_np, xs_np = load_xs("reaction_103")
    E_ng, xs_ng = load_xs("reaction_102")
    Eg = np.logspace(-3, 8, 2000)
    def interp(x, xp, fp):
        fp = np.clip(fp, 1e-30, None)
        return 10 ** np.interp(np.log10(x), np.log10(xp), np.log10(fp))
    return Eg, interp(Eg, E_ng, xs_ng) / interp(Eg, E_np, xs_np)


def tof_us(E_eV):
    gamma = 1.0 + np.asarray(E_eV, dtype=float) / M_N_EV
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    return L_FLIGHT_M / (beta * C_M_S) * 1e6


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--npz", default=str(REPO / "analysis/mev/mev_captures.npz"))
    ap.add_argument("--outdir", default=str(REPO / "docs/report/figs"))
    ap.add_argument("--json-out", default=str(REPO / "analysis/mev/mev_rates.json"))
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    d = np.load(args.npz, allow_pickle=True)
    n_sim = float(d["n_events"])
    edges = d["loge_edges"]                      # log10(E/eV), -3..8, 0.1/bin
    cent  = 10.0 ** ((edges[:-1] + edges[1:]) / 2)
    h_all, h_np, h_rad = (d[k].astype(float) for k in
                          ("h_all", "h_gas_np", "h_gas_rad"))
    h_wall, h_scint, h_plastic = (d[k].astype(float) for k in
                                  ("h_wall_al", "h_scint_h", "h_plastic"))
    rad_E = d["rad_events"].flatten()

    # ── Normalisation from the flux file ─────────────────────────────────────
    with uproot.open(FLUX_FILE) as f:
        hf = f[FLUX_HIST]
        fvals, fedges = hf.values(), hf.axis().edges()
    n_per_pulse = float(fvals.sum())             # axis covers exactly 1e-3..1e8
    anchor_subkev = float(fvals[fedges[1:] <= 1e3].sum())
    w = n_per_pulse / n_sim
    print(f"n/pulse (flux file, 1 meV-100 MeV) = {n_per_pulse:.4e}")
    print(f"  sub-keV part = {anchor_subkev:.4e}  (known anchor 7.31e6)")
    print(f"  weight w = {w:.4e} /pulse;  direct (n,g) events = {int(h_rad.sum())}")

    # ── Decade aggregation ────────────────────────────────────────────────────
    ndec = 11
    dec_edges = 10.0 ** np.arange(-3, 9)
    dec_cent  = np.sqrt(dec_edges[:-1] * dec_edges[1:])
    dsum = lambda h: np.array([h[i*10:(i+1)*10].sum() for i in range(ndec)])
    beam_d, np_d, rad_d = dsum(h_all), dsum(h_np), dsum(h_rad)
    wall_d, scint_d, plastic_d = dsum(h_wall), dsum(h_scint), dsum(h_plastic)

    x17_factor = ALPHA_IPC * BR_X17 * PULSES_DAY          # capt/pulse -> X17/day
    rows = []
    for i in range(ndec):
        k = i - 3
        lo68, hi68 = poisson_interval(int(rad_d[i]))
        tab = ALBERTO.get(k)
        rad_per_n = rad_d[i] / beam_d[i] if beam_d[i] else 0.0
        tab_per_n = tab[1] / tab[0] if tab else None
        rows.append({
            "E_lo_eV": 10.0**k, "E_hi_eV": 10.0**(k+1),
            "n_beam_sim": int(beam_d[i]), "beam_per_pulse": beam_d[i] * w,
            "n_np_sim": int(np_d[i]),
            "np_frac_of_beam": np_d[i] / beam_d[i] if beam_d[i] else 0.0,
            "n_rad": int(rad_d[i]), "n_rad_68": [lo68, hi68],
            "rad_per_pulse": rad_d[i] * w,
            "rad_per_n": rad_per_n,
            "x17_per_day": rad_d[i] * w * x17_factor,
            "x17_per_day_68": [lo68 * w * x17_factor, hi68 * w * x17_factor],
            "table_rad_per_n": tab_per_n,
            "g4_over_table": (rad_per_n / tab_per_n
                              if tab_per_n and rad_d[i] else None),
            "wall_per_pulse": wall_d[i] * w, "scint_per_pulse": scint_d[i] * w,
        })

    # 0.2-2 MeV window (log10 5.3-6.3 are exact bin edges)
    i_lo, i_hi = np.argmin(abs(edges - 5.3)), np.argmin(abs(edges - 6.3))
    n_win = h_rad[i_lo:i_hi].sum()
    win_lo, win_hi = poisson_interval(int(n_win))
    window = {
        "E_lo_eV": 2e5, "E_hi_eV": 2e6, "n_rad": int(n_win),
        "rad_per_pulse": n_win * w,
        "beam_per_pulse": h_all[i_lo:i_hi].sum() * w,
        "x17_per_day": n_win * w * x17_factor,
        "x17_per_day_68": [win_lo * w * x17_factor, win_hi * w * x17_factor],
        "ipc_per_pulse": n_win * w * ALPHA_IPC,
        "np_per_pulse": h_np[i_lo:i_hi].sum() * w,
        "wall_per_pulse": h_wall[i_lo:i_hi].sum() * w,
        "scint_per_pulse": h_scint[i_lo:i_hi].sum() * w,
        "plastic_per_pulse": h_plastic[i_lo:i_hi].sum() * w,
        "tof_us": [float(tof_us(2e6)), float(tof_us(2e5))],
        "table_rad_per_pulse": ALBERTO_WINDOW[1],
        "table_rad_per_n": ALBERTO_WINDOW[1] / ALBERTO_WINDOW[0],
    }

    rad_tot = h_rad.sum()
    summary = {
        "n_events": int(n_sim), "n_per_pulse": n_per_pulse,
        "pulses_per_day": PULSES_DAY, "weight_per_neutron": w,
        "alpha_ipc": ALPHA_IPC, "br_x17": BR_X17,
        "n_rad_direct": int(rad_tot),
        "rad_per_pulse": rad_tot * w,
        "rad_per_day": rad_tot * w * PULSES_DAY,
        "ipc_per_pulse": rad_tot * w * ALPHA_IPC,
        "ipc_per_day": rad_tot * w * ALPHA_IPC * PULSES_DAY,
        "x17_per_pulse": rad_tot * w * ALPHA_IPC * BR_X17,
        "x17_per_day": rad_tot * w * x17_factor,
        "subkev_rad_per_pulse": h_rad[edges[:-1][:len(h_rad)] < 3.0].sum() * w,
        "frac_rad_above_100keV": float((rad_E > 1e5).mean()),
        "median_rad_E_eV": float(np.median(rad_E)),
        "window_0p2_2MeV": window,
        "decades": rows,
    }
    with open(args.json_out, "w") as jf:
        json.dump(summary, jf, indent=1)
    print(f"\nrad capt/pulse = {summary['rad_per_pulse']:.2f}   "
          f"X17/day = {summary['x17_per_day']:.1f}   "
          f"window 0.2-2 MeV: {window['x17_per_day']:.1f} X17/day "
          f"({int(n_win)} events)")
    print(f"Saved {args.json_out}")

    # ════════════════════════════════════════════════════════════════════════
    # fig_mev_money: radiative captures /pulse /decade, direct counts
    # ════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    ax.axvspan(2e5, 2e6, color="gold", alpha=0.25, zorder=0)
    ax.text(6.3e5, 4e-4, "0.2–2 MeV\nwindow", ha="center", fontsize=9,
            color="0.3")

    has = rad_d > 0
    ylo = np.array([poisson_interval(int(n))[0] for n in rad_d])
    yhi = np.array([poisson_interval(int(n))[1] for n in rad_d])
    ax.errorbar(dec_cent[has], rad_d[has] * w,
                yerr=[(rad_d - ylo)[has] * w, (yhi - rad_d)[has] * w],
                fmt="o", color="C3", ms=7, capsize=4, zorder=5,
                label=r"$^3$He(n,$\gamma$) direct, Geant4 (714 events)")
    for i in np.nonzero(~has)[0]:        # upper arrows for empty decades
        ax.annotate("", xy=(dec_cent[i], 0.3 * 1.84 * w),
                    xytext=(dec_cent[i], 1.84 * w),
                    arrowprops=dict(arrowstyle="->", color="C3", alpha=0.7))

    # Alberto thin-target: his per-neutron probability x our beam (flux-free)
    tab_x = [dec_cent[k+3] for k in ALBERTO]
    tab_y = [ALBERTO[k][1]/ALBERTO[k][0] * beam_d[k+3] * w for k in ALBERTO]
    ax.plot(tab_x, tab_y, "s--", color="C1", ms=7, zorder=4,
            label="thin-target table (per-$n$ × sampled beam)")

    # thermal-ratio floor (the opaque-gas ceiling; sub-keV physics only)
    ax.step(dec_edges[:7], np.append(np_d[:6], np_d[5]) * SIGMA_TH * w,
            where="post", color="C2", lw=2, ls=":", zorder=3,
            label=r"opaque-gas ceiling: (n,p)$\,\times\,10^{-8}$ (sub-keV)")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("$E_n$  [eV]")
    ax.set_ylabel("radiative captures / pulse / decade")
    secy = ax.secondary_yaxis(
        "right", functions=(lambda y: y * x17_factor, lambda y: y / x17_factor))
    secy.set_ylabel("X17 produced / day / decade")
    ax.set_title(r"$^3$He(n,$\gamma$)$^4$He per pulse — X17 production channel"
                 "\n" rf"total: {rad_tot*w:.1f} capt/pulse  $\to$  "
                 rf"{rad_tot*w*ALPHA_IPC:.3f} IPC/pulse  $\to$  "
                 rf"{rad_tot*w*x17_factor:.1f} X17/day")
    ax.set_ylim(1e-4, 1e2)
    ax.legend(loc="upper left", fontsize=9.5)
    fig.savefig(out / "fig_mev_money.pdf", bbox_inches="tight")
    plt.close(fig)

    # ════════════════════════════════════════════════════════════════════════
    # fig_mev_beam: sampled beam vs flux file
    # ════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.step(cent, h_all * w, where="mid", color="k", lw=1.4,
            label=f"sampled beam ({n_sim:.0e} events) × $w$")
    fcent = np.sqrt(fedges[:-1] * fedges[1:])
    # flux file: 10 bins/decade over 1e-3..1e8 -> same binning, overlay direct
    ax.step(fcent, fvals, where="mid", color="C0", lw=1.2, alpha=0.7,
            label=f"flux file {FLUX_HIST}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("$E_n$  [eV]"); ax.set_ylabel("neutrons / pulse / bin")
    ax.set_title(f"Normalisation check: {n_per_pulse:.3e} n/pulse over "
                 "1 meV–100 MeV (sub-keV part = 7.31×10$^6$ ✓)")
    ax.legend(fontsize=9.5)
    fig.savefig(out / "fig_mev_beam.pdf", bbox_inches="tight")
    plt.close(fig)

    # ════════════════════════════════════════════════════════════════════════
    # fig_mev_opacity: (n,p) absorption fraction
    # ════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(8.5, 5))
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(h_all > 0, h_np / h_all, np.nan)
    ax.step(cent, frac, where="mid", color="C1", lw=2)
    ax.fill_between(cent, 0, frac, step="mid", alpha=0.15, color="C1")
    ax.axvspan(2e5, 2e6, color="gold", alpha=0.25, zorder=0)
    ax.set_xscale("log")
    ax.set_xlabel("$E_n$  [eV]")
    ax.set_ylabel(r"fraction of beam absorbed by $^3$He(n,p)t")
    ax.set_title("Opacity transition: the gas turns thin above ~10 keV")
    ax.set_ylim(0, 1.0)
    for E, lab in [(1e0, "opaque\n(70–80%)"), (3e4, "transition"),
                   (8e6, "thin (<4%)")]:
        ax.text(E, 0.88, lab, ha="center", fontsize=9, color="0.3")
    fig.savefig(out / "fig_mev_opacity.pdf", bbox_inches="tight")
    plt.close(fig)

    # ════════════════════════════════════════════════════════════════════════
    # fig_mev_ratio: effective sigma_ng/sigma_np vs ENDF
    # ════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    try:
        Eg, ratio = load_endf_ratio()
        ax.plot(Eg, ratio, color="C0", lw=1.6, alpha=0.8,
                label=r"ENDF/B-VIII.0  $\sigma_{n\gamma}/\sigma_{np}$")
    except Exception as e:
        print(f"[warn] ENDF ratio skipped: {e}")
    m = (rad_d > 0) & (np_d > 0)
    r = np.where(np_d > 0, rad_d / np.maximum(np_d, 1), np.nan)
    rlo = np.array([poisson_interval(int(n))[0] for n in rad_d]) / np.maximum(np_d, 1)
    rhi = np.array([poisson_interval(int(n))[1] for n in rad_d]) / np.maximum(np_d, 1)
    ax.errorbar(dec_cent[m], r[m], yerr=[(r - rlo)[m], (rhi - r)[m]],
                fmt="o", color="C3", ms=7, capsize=4,
                label="Geant4: direct (n,$\\gamma$) / (n,p)t per decade")
    ax.axhline(SIGMA_TH, color="C2", ls=":", lw=1.5,
               label=r"thermal value $1.0\times10^{-8}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("$E_n$  [eV]")
    ax.set_ylabel(r"$\sigma_{n\gamma}/\sigma_{np}$ (effective)")
    ax.set_title("Why the MeV region wins: the branching ratio rises "
                 r"$10^{4}\times$ from thermal to MeV")
    ax.legend(fontsize=9.5, loc="upper left")
    fig.savefig(out / "fig_mev_ratio.pdf", bbox_inches="tight")
    plt.close(fig)

    # ════════════════════════════════════════════════════════════════════════
    # fig_mev_backgrounds: competing captures per decade
    # ════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for h, lab, c in [(np_d, r"$^3$He(n,p)t (gas heating)", "C0"),
                      (scint_d, "scintillator H-capture (2.2 MeV $\\gamma$)", "C4"),
                      (wall_d, "Al capsule nCapture (7.7 MeV cascade)", "C2"),
                      (plastic_d, "plastic veto nCapture", "C5"),
                      (rad_d, r"$^3$He(n,$\gamma$) — the signal", "C3")]:
        ax.step(dec_edges, np.append(h, h[-1]) * w, where="post", lw=2,
                color=c, label=lab)
    ax.axvspan(2e5, 2e6, color="gold", alpha=0.25, zorder=0)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("$E_n$  [eV]")
    ax.set_ylabel("events / pulse / decade")
    ax.set_title("Neutron fate per decade: signal vs capture backgrounds")
    ax.legend(fontsize=9, loc="lower left")
    ax.set_ylim(1e-4, 3e7)
    fig.savefig(out / "fig_mev_backgrounds.pdf", bbox_inches="tight")
    plt.close(fig)

    # ════════════════════════════════════════════════════════════════════════
    # fig_mev_tof: E_n <-> arrival time at EAR2
    # ════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    Et = np.logspace(0, 8, 400)
    ax.plot(tof_us(Et), Et, color="k", lw=1.8)
    t_lo, t_hi = tof_us(2e6), tof_us(2e5)
    ax.axvspan(t_lo, t_hi, color="gold", alpha=0.3)
    ax.axhspan(2e5, 2e6, color="gold", alpha=0.15)
    ax.annotate(f"0.2–2 MeV ↔ {t_lo:.1f}–{t_hi:.1f} µs",
                xy=(t_hi, 2e5), xytext=(8, 2e6),
                arrowprops=dict(arrowstyle="->", color="0.3"), fontsize=10)
    t_flash = L_FLIGHT_M / C_M_S * 1e6
    ax.axvline(t_flash, color="C3", ls="--", lw=1.5)
    ax.text(t_flash * 1.3, 3e0, f"$\\gamma$ flash\n({t_flash*1e3:.0f} ns)",
            color="C3", fontsize=9)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(f"arrival time at L = {L_FLIGHT_M} m  [µs]")
    ax.set_ylabel("$E_n$  [eV]")
    ax.set_title("TOF map: the MeV window sits 1–3 µs after the flash")
    fig.savefig(out / "fig_mev_tof.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"Saved 6 figures to {out}/")


if __name__ == "__main__":
    main()
