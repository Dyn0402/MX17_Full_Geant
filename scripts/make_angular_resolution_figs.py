#!/usr/bin/env python3
"""make_angular_resolution_figs.py — figures for docs/angular_resolution/

Reads analysis/pairs_v2/geant4_response.json (10M-event pairs_v2 run, STEP
small-capsule geometry) and produces the four figures referenced by
docs/angular_resolution/angular_resolution_note.md:

  fig_psi_vs_ke.png        median per-particle direction error vs KE,
                           per estimator, + Highland MS prediction
  fig_sigma68_vs_minke.png pair opening-angle sigma68 + bias vs min(KE-,KE+)
  fig_ms_budget.png        upstream material budget (new STEP geometry)
  fig_theta_dilution.png   truth at-rest X17 opening-angle spectrum vs the
                           same spectrum smeared by the measured resolution

Usage:  python scripts/make_angular_resolution_figs.py
"""

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
RESP = REPO / "analysis/pairs_v2/geant4_response.json"
OUT  = REPO / "docs/angular_resolution/figs"

ME = 0.511          # MeV
M_X17 = 16.8        # MeV
E_X17 = 20.58       # MeV total energy of at-rest-capture X17

# ── Upstream material budget, STEP capsule geometry (current sim) ────────────
# (label, thickness [cm] at normal incidence, X0 [cm])
# He-3 X0: 94.32 g/cm2 (He-4) scaled by A: 3/4 -> 70.7 g/cm2; rho = 62.7 mg/cm3.
MS_LAYERS = [
    ("He-3 gas 500 bar (~10path)", 1.0,    70.7 / 62.7e-3),
    ("Al barrel wall 0.6 mm",          0.06,   8.897),
    ("CFRP wrap 0.9 mm",               0.09,   42.70 / 1.55),
    ("Air gap ~23.9 cm",               23.85,  36.62 / 1.205e-3),
    ("Mylar window 40 um",             0.004,  28.54),
    ("Kapton cathode 50 um",           0.005,  28.58),
    ("Cu cathode 9 um",                0.0009, 1.436),
]
COLORS = ["#4a90d9", "#d62728", "#ff7f0e", "#2ca02c",
          "#9467bd", "#8c564b", "#e377c2"]


def highland(p_mev, x_over_x0):
    """Plane-projected theta0 [rad] for electrons of momentum p."""
    x = max(x_over_x0, 1e-12)
    beta = p_mev / np.sqrt(p_mev**2 + ME**2)
    return (13.6 / (beta * p_mev)) * np.sqrt(x) * (1 + 0.038 * np.log(x))


def hist_quantile(h, cen, q):
    n = h.sum()
    if n < 30:
        return np.nan
    cdf = np.cumsum(h) / n
    return cen[min(np.searchsorted(cdf, q), len(cen) - 1)]


def fig_psi_vs_ke(d, out=None, ext="png"):
    out = out or OUT
    de = d["direction_error"]
    ke_e  = np.array(de["ke_bin_edges"]);  ke_c  = 0.5 * (ke_e[:-1] + ke_e[1:])
    psi_e = np.array(de["psi_bin_edges"]); psi_c = 0.5 * (psi_e[:-1] + psi_e[1:])

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    style = {"first":   ("#000000", "-",  "dir @ 1st MM hit (perfect tracker)"),
             "fit":     ("#2ca02c", "--", "PCA fit of drift hits"),
             "vline":   ("#1f77b4", "-",  "true vtx → 1st MM hit (ideal chord)"),
             "nomline": ("#9467bd", "-.", "target centre → 1st MM hit (realistic chord)")}
    for nm, (col, ls, lbl) in style.items():
        h = np.array(de["estimators"][nm]["counts"], dtype=float)
        med = np.array([hist_quantile(row, psi_c, 0.5) for row in h])
        v = ~np.isnan(med)
        ax.plot(ke_c[v], med[v], color=col, ls=ls, lw=2, label=lbl)

    tot = sum(t / x0 for _, t, x0 in MS_LAYERS)
    kes = np.linspace(0.7, 20, 200)
    ps  = np.sqrt((kes + ME)**2 - ME**2)
    ax.plot(kes, 1.177 * np.degrees([highland(p, tot) for p in ps]),
            "r:", lw=2.2,
            label=f"Highland upstream total (x/X0={tot*100:.2f}%), median space angle")

    ax.axhline(45, color="grey", lw=0.8, ls=":")
    ax.text(15, 43, "histogram cap (45°)", fontsize=8, color="grey", va="top")
    ax.set_xlabel("True particle KE  [MeV]")
    ax.set_ylabel("Median ψ vs truth direction at vertex  [deg]")
    ax.set_ylim(0, 47); ax.set_xlim(0, 20)
    ax.set_title("Per-particle direction error vs KE — pairs_v2 (10M events, "
                 "STEP capsule)\nall estimators converge to the upstream-MS floor")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(out / f"fig_psi_vs_ke.{ext}", dpi=150)
    plt.close(fig)


def fig_sigma68_vs_minke(d, out=None, ext="png"):
    out = out or OUT
    v = d["validation"]
    edges = np.array(v["minke_bin_edges"]); cen = 0.5 * (edges[:-1] + edges[1:])
    style = {"first": ("#000000", "dir @ 1st MM hit"),
             "vline": ("#1f77b4", "true vtx → 1st MM hit (ideal)"),
             "fit":   ("#2ca02c", "PCA drift-track fit")}

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for mth, (col, lbl) in style.items():
        s68  = np.array([np.nan if x is None else x
                         for x in v["methods"][mth]["sigma68_deg"]])
        bias = np.array([np.nan if x is None else x
                         for x in v["methods"][mth]["bias_deg"]])
        ok = ~np.isnan(s68)
        axes[0].plot(cen[ok], s68[ok],  color=col, lw=2, marker="o", ms=3, label=lbl)
        axes[1].plot(cen[ok], bias[ok], color=col, lw=2, marker="o", ms=3, label=lbl)
    axes[0].set_ylabel("σ68(Δθ)  [deg]"); axes[0].set_ylim(0, 45)
    axes[0].set_title("Opening-angle resolution vs pair symmetry")
    axes[1].set_ylabel("median Δθ (bias)  [deg]")
    axes[1].axhline(0, color="grey", lw=0.8, ls="--")
    axes[1].set_title("Bias (correctable on average)")
    for ax in axes:
        ax.set_xlabel("min(KE⁻, KE⁺)  [MeV]")
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3); ax.set_xlim(0, 11)
    fig.suptitle("Pair opening-angle Δθ = θ_reco − θ_truth vs min pair KE  "
                 "(X17 + IPC combined)", fontsize=12)
    fig.tight_layout(); fig.savefig(out / f"fig_sigma68_vs_minke.{ext}", dpi=150)
    plt.close(fig)


def fig_ms_budget(out=None, ext="png"):
    out = out or OUT
    names = [l[0] for l in MS_LAYERS]
    xfrac = [t / x0 for _, t, x0 in MS_LAYERS]
    tot = sum(xfrac)
    p_ref = 10.0

    fig, ax = plt.subplots(figsize=(9, 4.8))
    y = np.arange(len(names))
    ax.barh(y, np.array(xfrac) * 100, color=COLORS, edgecolor="k", lw=0.4)
    for yi, xf in zip(y, xfrac):
        ax.text(xf * 100 * 1.12, yi,
                f"θ₀={np.degrees(highland(p_ref, xf)):.2f}°   ({xf/tot*100:.0f}%)",
                va="center", fontsize=8.5)
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=9)
    ax.invert_yaxis(); ax.set_xscale("log"); ax.set_xlim(right=10)
    ax.set_xlabel("x / X₀  [%]  (normal incidence)")
    ax.set_title(f"Upstream material budget, STEP capsule geometry — "
                 f"total x/X₀ = {tot*100:.2f}%\n"
                 f"θ₀ = {np.degrees(highland(p_ref, tot)):.1f}° per particle at "
                 f"p = {p_ref:.0f} MeV/c  (capsule wall = "
                 f"{(xfrac[1]+xfrac[2])/tot*100:.0f}% of budget)")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout(); fig.savefig(out / f"fig_ms_budget.{ext}", dpi=150)
    plt.close(fig)


def sample_pairs(m_parent, n, rng):
    """Full lab pair kinematics for e+e- pairs from a parent of mass m_parent
    (scalar or per-event array) carrying total energy E* = 20.58 MeV, emitted
    isotropically from an at-rest 4He* and decaying isotropically — the exact
    kinematics of GeneratePairFromParent in X17PrimaryGenerator.cc.

    Returns (theta_deg, ke1, ke2, u1, u2): opening angle, per-leg lab kinetic
    energies [MeV], per-leg unit direction vectors (n,3)."""
    m = np.broadcast_to(np.asarray(m_parent, dtype=float), (n,))
    gam, betgam = E_X17 / m, np.sqrt(E_X17**2 - m**2) / m
    e_cm = m / 2
    p_cm = np.sqrt(np.maximum(e_cm**2 - ME**2, 0.0))
    cth = rng.uniform(-1, 1, n); sth = np.sqrt(1 - cth**2)
    phi = rng.uniform(0, 2 * np.pi, n)
    # lepton momenta in parent frame (z = boost axis); e+ is -p
    pz, px, py = p_cm * cth, p_cm * sth * np.cos(phi), p_cm * sth * np.sin(phi)
    pz1 = gam * pz + betgam * e_cm;  pz2 = -gam * pz + betgam * e_cm
    ke1 = gam * e_cm + betgam * pz - ME
    ke2 = gam * e_cm - betgam * pz - ME
    p1 = np.stack([px, py, pz1], 1); p2 = np.stack([-px, -py, pz2], 1)
    u1 = p1 / np.linalg.norm(p1, axis=1, keepdims=True)
    u2 = p2 / np.linalg.norm(p2, axis=1, keepdims=True)
    cos_open = (u1 * u2).sum(1)
    th = np.degrees(np.arccos(np.clip(cos_open, -1, 1)))
    return th, ke1, ke2, u1, u2


def sample_pair_angles(m_parent, n, rng):
    """Lab opening angles [deg] only (see sample_pairs)."""
    return sample_pairs(m_parent, n, rng)[0]


def sample_ipc_masses(n, rng):
    """IPC virtual-state masses: dN/dM ∝ 1/M on [2me, E*] (generator's
    inverse-CDF sampling)."""
    return 2 * ME * (E_X17 / (2 * ME)) ** rng.uniform(0, 1, n)


def smear_angles(th, s68, bias, rng):
    """Gaussian opening-angle smear with reflection at the physical
    boundaries 0° and 180° (an opening angle cannot leave [0, 180])."""
    t = th + rng.normal(bias, s68, len(th))
    t = np.abs(t)
    return np.where(t > 180.0, 360.0 - t, t)


def make_psi_sampler(d, estimator="nomline"):
    """Per-leg direction-error sampler from the measured Geant4 response.

    Returns sample(ke, rng) → ψ [rad] drawn from P(ψ | KE) of the given
    estimator in geant4_response.json.  Draws landing in the 44.75–45° cap
    bin (overflow: ψ exceeded the histogram range, soft legs) are spread
    uniformly over 45–90° — a broad-tail approximation."""
    de = d["direction_error"]
    ke_edges  = np.array(de["ke_bin_edges"])
    psi_edges = np.array(de["psi_bin_edges"])
    counts = np.array(de["estimators"][estimator]["counts"], float)
    row_n = counts.sum(1)
    valid = np.where(row_n > 0)[0]
    for i in np.where(row_n == 0)[0]:          # empty KE rows → nearest
        counts[i] = counts[valid[np.argmin(np.abs(valid - i))]]
    cdf = np.cumsum(counts, axis=1)
    cdf /= cdf[:, -1:]
    n_psi = len(psi_edges) - 1
    dpsi = np.diff(psi_edges)

    def sample(ke, rng):
        ki = np.clip(np.digitize(ke, ke_edges) - 1, 0, len(ke_edges) - 2)
        u = rng.uniform(0, 1, len(ke))
        pi = np.clip((cdf[ki] < u[:, None]).sum(1), 0, n_psi - 1)
        psi = psi_edges[pi] + rng.uniform(0, 1, len(ke)) * dpsi[pi]
        cap = pi == n_psi - 1                  # overflow bin → 45–90° tail
        psi[cap] = rng.uniform(45.0, 90.0, cap.sum())
        return np.radians(psi)

    return sample


def smear_pair_response(ke1, ke2, u1, u2, sampler, rng):
    """Opening angles [deg] after smearing each leg's direction by a ψ drawn
    from the measured P(ψ | KE), with uniform azimuth around the true
    direction — reproduces the pair Δθ resolution, its θ/energy dependence,
    and the positive bias, per the geant4_response.json design notes."""
    def tilt(u, ke):
        psi = sampler(ke, rng)
        a = np.zeros_like(u); a[:, 2] = 1.0
        near_z = np.abs(u[:, 2]) > 0.9
        a[near_z] = np.array([1.0, 0.0, 0.0])
        e1 = np.cross(u, a); e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
        e2 = np.cross(u, e1)
        phi = rng.uniform(0, 2 * np.pi, len(u))
        return (np.cos(psi)[:, None] * u +
                np.sin(psi)[:, None] * (np.cos(phi)[:, None] * e1 +
                                        np.sin(phi)[:, None] * e2))
    c = (tilt(u1, ke1) * tilt(u2, ke2)).sum(1)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def fig_theta_dilution(d, n=2_000_000, seed=42, out=None, ext="png"):
    """Truth at-rest X17 opening angle vs the same smeared by the measured
    Δθ response (Gaussian, KE-independent approximation: σ = σ68, +bias)."""
    out = out or OUT
    rng = np.random.default_rng(seed)
    th = sample_pair_angles(M_X17, n, rng)

    bins = np.linspace(0, 180, 91)
    cen = 0.5 * (bins[:-1] + bins[1:])
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    h0, _ = np.histogram(th, bins=bins)
    ax.step(cen, h0 / h0.sum(), where="mid", color="#2ca02c", lw=2,
            label="truth (at-rest capture, E*=20.58 MeV)")
    for s68, bias, col, lbl in [
        (15.5, 3.0, "#e84040", "smeared σ68=15.5°, bias +3° (full spectrum)"),
        (12.0, 2.5, "#1f77b4", "smeared σ68=12°, bias +2.5° (symmetric pairs, ideal chord)"),
    ]:
        hs, _ = np.histogram(smear_angles(th, s68, bias, rng), bins=bins)
        ax.step(cen, hs / hs.sum(), where="mid", color=col, lw=2, ls="--",
                label=lbl)
    ax.set_xlabel("e⁺e⁻ opening angle  [deg]")
    ax.set_ylabel("normalised counts / 2°")
    ax.set_title("X17 opening-angle peak dilution by the measured resolution\n"
                 "(acceptance not applied; Gaussian-response illustration)")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3); ax.set_xlim(60, 180)
    fig.tight_layout(); fig.savefig(out / f"fig_theta_dilution.{ext}", dpi=150)
    plt.close(fig)


# Recorded yields for a 30-day run, ~10 us per-pulse readout (no scint
# trigger, 2026-06-12 decision): the readout covers E_n >~ 20 keV, i.e.
# essentially full production (32.4 X17/day, mev_rates.json decades >=10 keV),
# and every pair in MM geometric acceptance is recorded ("MM double",
# pairs_v2: 19.6% X17 / 23.6% IPC — same-arm pairs now count).
# X17: 32.4/day x 19.6% = 6.4/day; IPC: (32.4/0.025)/day x 23.6% = 306/day.
N_X17_30D = 191.0
N_IPC_30D = 9190.0
# Gaussian-response smear parameters (sigma68, bias) per species.
# SMEAR_* : current default method (dir @ 1st MM hit), full spectrum.
# BEST_*  : best feasible — target-centre chord ("nomline": project back to
#           the capsule centre), full spectrum.  Used for the money plot.
SMEAR_X17 = (15.5, 3.0)
SMEAR_IPC = (15.0, 3.0)
BEST_X17  = (15.0, 2.5)
BEST_IPC  = (13.5, 2.5)


def fig_ipc_dilution(n=2_000_000, seed=43, out=None, ext="png"):
    """Truth IPC opening-angle spectrum vs the same smeared by the measured
    resolution — companion to fig_theta_dilution."""
    out = out or OUT
    rng = np.random.default_rng(seed)
    th = sample_pair_angles(sample_ipc_masses(n, rng), n, rng)

    bins = np.linspace(0, 180, 91)
    cen = 0.5 * (bins[:-1] + bins[1:])
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    h0, _ = np.histogram(th, bins=bins)
    ax.step(cen, h0 / h0.sum(), where="mid", color="#4a90d9", lw=2,
            label="truth (IPC, dN/dM ∝ 1/M, E*=20.58 MeV)")
    for s68, bias, col, lbl in [
        (SMEAR_IPC[0], SMEAR_IPC[1], "#e84040",
         f"smeared σ68={SMEAR_IPC[0]:.0f}°, bias +{SMEAR_IPC[1]:.0f}° (full spectrum)"),
        (13.0, 2.5, "#1f77b4", "smeared σ68=13°, bias +2.5° (ideal chord)"),
    ]:
        hs, _ = np.histogram(smear_angles(th, s68, bias, rng), bins=bins)
        ax.step(cen, hs / hs.sum(), where="mid", color=col, lw=2, ls="--",
                label=lbl)
    ax.set_xlabel("e⁺e⁻ opening angle  [deg]")
    ax.set_ylabel("normalised counts / 2°")
    ax.set_title("IPC background: opening-angle spectrum before and after the "
                 "measured resolution\n(acceptance not applied; "
                 "Gaussian-response illustration)")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3); ax.set_xlim(0, 180)
    fig.tight_layout(); fig.savefig(out / f"fig_ipc_dilution.{ext}", dpi=150)
    plt.close(fig)


def fig_theta_money(d, n=2_000_000, seed=44, out=None, ext="png"):
    """The money plot: the truth distribution vs the best we can do given
    capsule multiple scattering.

    Left   : X17 signal alone — truth opening-angle spectrum vs the same
             after per-leg response smearing with the best feasible
             estimator (target-centre chord, "nomline").
    Middle : X17 stacked on the IPC continuum, truth shapes.
    Right  : the same, both smeared — the spectrum the analysis actually
             has to extract the X17 signal from.

    Smearing: each leg's direction error psi is drawn from the measured
    Geant4 P(psi | KE) tables, so the pair resolution varies with opening
    angle and species (best at the X17 shoulder where pairs are symmetric).

    Normalisation: 30-day recorded yields under the 10 us readout +
    MM-double acceptance (191 X17 / 9190 IPC).  This is a BEST-CASE picture:
    smallest feasible smearing, no MM acceptance shaping vs theta, no
    pile-up backgrounds, no gamma-flash losses."""
    out = out or OUT
    rng = np.random.default_rng(seed)
    sampler = make_psi_sampler(d, "nomline")
    th_x, ke1x, ke2x, u1x, u2x = sample_pairs(M_X17, n, rng)
    th_i, ke1i, ke2i, u1i, u2i = sample_pairs(sample_ipc_masses(n, rng), n, rng)
    th_x_s = smear_pair_response(ke1x, ke2x, u1x, u2x, sampler, rng)
    th_i_s = smear_pair_response(ke1i, ke2i, u1i, u2i, sampler, rng)

    bins = np.linspace(0, 180, 91)
    cen = 0.5 * (bins[:-1] + bins[1:])

    def counts(th, total):
        h, _ = np.histogram(th, bins=bins)
        return h / h.sum() * total

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4))

    # ── left: the signal alone, expected vs observable ───────────────────
    ax = axes[0]
    sx_t, sx_s = counts(th_x, N_X17_30D), counts(th_x_s, N_X17_30D)
    ax.fill_between(cen, 0, sx_t, step="mid", color="#2ca02c", alpha=0.25)
    ax.step(cen, sx_t, where="mid", color="#2ca02c", lw=2,
            label="truth distribution")
    ax.step(cen, sx_s, where="mid", color="#e84040", lw=2, ls="--",
            label="measured response — best we can get")
    ax.set_ylabel("expected events / 2°  (30 days, 10 µs readout, "
                  "MM acceptance)")
    ax.set_title(f"X17 signal alone  ({N_X17_30D:.0f} events)")

    # ── middle/right: stacked on IPC, truth then smeared ─────────────────
    for ax, b, s, title in [
        (axes[1], counts(th_i, N_IPC_30D), sx_t,
         "stacked on IPC — truth shapes"),
        (axes[2], counts(th_i_s, N_IPC_30D), sx_s,
         "stacked + smeared — what we must analyse"),
    ]:
        ax.fill_between(cen, 0, b, step="mid", color="#9bb8d4", alpha=0.9,
                        label=f"IPC  ({N_IPC_30D:.0f} events)")
        ax.fill_between(cen, b, b + s, step="mid", color="#e84040", alpha=0.9,
                        label=f"X17 stacked  ({N_X17_30D:.0f} events)")
        ax.step(cen, b + s, where="mid", color="k", lw=1.0)
        ax.set_title(title)

    axes[2].step(cen, sx_s, where="mid", color="#e84040", lw=2, ls="--",
            label=f"X17 ({N_X17_30D:.0f} events)")

    for ax in axes:
        ax.set_xlabel("e⁺e⁻ opening angle  [deg]")
        loc = "upper left" if ax is axes[0] else "upper right"
        ax.legend(fontsize=9, loc=loc); ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 180); ax.set_ylim(bottom=0)

    fig.suptitle("Best we can do with capsule multiple scattering — "
                 "30-day recorded yields (10 µs readout, MM-double "
                 "acceptance)\nBEST CASE: per-leg smearing from the measured "
                 "P(ψ|KE) response, best estimator (target-centre chord); "
                 "production-level shapes, no pile-up backgrounds, no "
                 "γ-flash losses", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / f"fig_theta_money.{ext}", dpi=150)
    plt.close(fig)


def fig_res_vs_theta(d, n=2_000_000, seed=45, out=None, ext="png"):
    """Pair opening-angle resolution and bias vs TRUTH opening angle, per
    species, from per-leg response smearing (best estimator).  Shows that the
    resolution is not flat: X17 is sharpest at the 109° shoulder (symmetric
    pairs) and degrades toward 180° (asymmetric, soft leg); IPC is sharpest
    at small angle (low-mass pairs split E* evenly → two stiff legs)."""
    out = out or OUT
    rng = np.random.default_rng(seed)
    sampler = make_psi_sampler(d, "nomline")
    edges = np.arange(0, 181, 7.5)
    cen = 0.5 * (edges[:-1] + edges[1:])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for lab, mass, col in [("X17", M_X17, "#e84040"),
                           ("IPC", None, "#4a90d9")]:
        m = mass if mass is not None else sample_ipc_masses(n, rng)
        th, ke1, ke2, u1, u2 = sample_pairs(m, n, rng)
        dth = smear_pair_response(ke1, ke2, u1, u2, sampler, rng) - th
        s68 = np.full(len(cen), np.nan); bias = np.full(len(cen), np.nan)
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            sel = (th > lo) & (th < hi)
            if sel.sum() < 3000:
                continue
            q16, q50, q84 = np.percentile(dth[sel], [16, 50, 84])
            s68[i] = 0.5 * (q84 - q16); bias[i] = q50
        ok = ~np.isnan(s68)
        axes[0].plot(cen[ok], s68[ok], color=col, lw=2, marker="o", ms=3,
                     label=lab)
        axes[1].plot(cen[ok], bias[ok], color=col, lw=2, marker="o", ms=3,
                     label=lab)
    for ax in axes:
        ax.axvline(109.4, color="#2ca02c", lw=1.2, ls=":",
                   label="X17 shoulder (109°)")
        ax.set_xlabel("truth e⁺e⁻ opening angle  [deg]")
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3); ax.set_xlim(0, 180)
    axes[0].set_ylabel("σ68(Δθ)  [deg]"); axes[0].set_ylim(bottom=0)
    axes[0].set_title("resolution vs truth angle")
    axes[1].set_ylabel("median Δθ (bias)  [deg]")
    axes[1].axhline(0, color="grey", lw=0.8, ls="--")
    axes[1].set_title("bias vs truth angle")
    fig.suptitle("Opening-angle response vs truth angle, per species — "
                 "per-leg P(ψ|KE) smearing, target-centre chord\n"
                 "(X17 sharpest at its shoulder: symmetric pairs; "
                 "boundary effects drive the bias at 0°/180°)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / f"fig_res_vs_theta.{ext}", dpi=150)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = json.load(open(RESP))
    fig_psi_vs_ke(d)
    fig_sigma68_vs_minke(d)
    fig_ms_budget()
    fig_theta_dilution(d)
    fig_ipc_dilution()
    fig_res_vs_theta(d)
    fig_theta_money(d)
    print(f"7 figures written to {OUT}")


if __name__ == "__main__":
    main()
