#!/usr/bin/env python3
"""
study_invmass.py — invariant-mass reconstruction strategy + resolution study
=============================================================================
Works out how to implement e+e- invariant-mass reconstruction with the LS
calorimeter stack and what resolution is realistically achievable.

Reads the same mx17_full_sim ROOT files as analyze_pairs.py, but caches a
per-event table (pickle) so the reconstruction study can iterate quickly.

Key physics exploited
---------------------
1. M² = 2mₑ² + 2(E₁E₂ − p₁p₂cosθ)  — needs two energies and the opening angle.
2. The ⁴He* transition fixes  E₁+E₂ = 20.58 MeV  (recoil ~keV) for BOTH X17
   and IPC events.  Writing E₁ = xS, E₂ = (1−x)S with S known, the calorimeter
   only has to measure the energy-sharing fraction x.  Since
   ∂M/∂x ∝ (1−2x) = 0 at x = ½, the mass is first-order INSENSITIVE to
   calorimeter errors for symmetric pairs (X17 peaks at x ≈ ½).
3. Experimentally we measure energy PER ARM, not per particle:
   the arm-sum estimator (all scintillator edep in the arm the MM track points
   to) is compared against the idealised ancestry-attributed sum.

Reconstruction ladder (X17 signal, diff-arm MM topology):
   truth E + truth θ        → sanity (= 16.8)
   truth E + reco θ         → angular term alone
   reco  E + truth θ        → calorimeter term alone
   reco  E + reco θ         → direct full reconstruction
   constraint (x, S) + reco θ → constrained full reconstruction
plus photostatistics scenarios σ_E/E = k/√E on top of the deposited energy.

Usage:
    python study_invmass.py ~/Desktop/x17_pairs/*.root -o invmass_study.pdf
    python study_invmass.py --cache-only ~/Desktop/x17_pairs/*.root
"""

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_pairs as ap

ME    = 0.511
E_SUM = 20.58          # MeV — total pair energy fixed by the 4He* transition
M_X17 = 16.8

SCINT_ALL = frozenset(["PlasticScint", "LiqScint_1", "LiqScint_2",
                       "BackScintL", "BackScintR"])

KEEP_COLS = (
    ["eventID", "event_type", "inv_mass", "em_ke", "ep_ke",
     "em_px", "em_py", "em_pz", "ep_px", "ep_py", "ep_pz", "openingAngle"]
    + [f"{p}_{c}" for p in ("em", "ep") for c in
       ["edep_ls", "edep_bs", "edep_all", "mm_arm",
        "mm_gx", "mm_gy", "mm_gz", "mm_px", "mm_py", "mm_pz",
        "trk_n", "trk_s_px", "trk_s_py", "trk_s_pz", "trig", "stop_ke"]]
    + ["single_trig", "double_trig", "mm_double", "mm_same_arm", "mm_diff_arm"]
    + [f"arm_{t}_{a}" for t in ("ls", "all") for a in range(4)]
)


# ── Per-event table construction ──────────────────────────────────────────────

def _arm_sums(df):
    """Per-event, per-arm scintillator edep sums [MeV] — the experimentally
    measurable calorimeter quantity (no ancestry attribution).
      arm_ls_<a>  : LiqScint_1+2 in arm a
      arm_all_<a> : PlasticScint + LS + BackScint in arm a
    """
    det  = ap._decode_col(df["detType"])
    eids = pd.Index(df["eventID"].unique(), name="eventID")
    out  = pd.DataFrame(index=eids)
    for tag, layers in [("ls", ap.LS_LAYERS), ("all", SCINT_ALL)]:
        sel = det.isin(layers)
        if sel.any():
            piv = (df.loc[sel].groupby(["eventID", "armID"])["edep"].sum()
                   / 1e6).unstack("armID")
        else:
            piv = pd.DataFrame(index=eids)
        for a in range(4):
            col = piv[a] if a in piv.columns else pd.Series(np.nan, index=piv.index)
            out[f"arm_{tag}_{a}"] = col.reindex(eids).fillna(0.0)
    return out.reset_index()


def process_file_events(filepath, chunk_size=300_000):
    """Like ap.process_file but returns the merged per-event DataFrame."""
    with uproot.open(filepath) as f:
        evt_df = f["EventTree"].arrays(library="pd")
        if evt_df.empty:
            return evt_df

        leftover, summaries, armsums = None, [], []
        for chunk in f["HitTree"].iterate(ap.HIT_COLS, step_size=chunk_size,
                                          library="pd"):
            chunk["detType"]  = np.asarray(chunk["detType"])
            chunk["particle"] = np.asarray(chunk["particle"])
            if leftover is not None:
                chunk = pd.concat([leftover, chunk], ignore_index=True)
            last_eid = chunk["eventID"].iloc[-1]
            complete = chunk[chunk["eventID"] < last_eid]
            leftover = chunk[chunk["eventID"] >= last_eid]
            if not complete.empty:
                summaries.append(ap._process_hits(complete.copy()))
                armsums.append(_arm_sums(complete))
            del chunk
        if leftover is not None and not leftover.empty:
            summaries.append(ap._process_hits(leftover.copy()))
            armsums.append(_arm_sums(leftover))

        if not summaries:
            return evt_df
        merged = (evt_df
                  .merge(pd.concat(summaries, ignore_index=True),
                         on="eventID", how="left")
                  .merge(pd.concat(armsums, ignore_index=True),
                         on="eventID", how="left"))
        return merged[[c for c in KEEP_COLS if c in merged.columns]]


def build_table(files, cache_path, rebuild=False, workers=4):
    if cache_path and os.path.exists(cache_path) and not rebuild:
        print(f"Loading cached per-event table: {cache_path}")
        return pd.read_pickle(cache_path)
    print(f"Processing {len(files)} files with {workers} workers ...")
    parts = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(process_file_events, fp): fp for fp in files}
        for fut in as_completed(futs):
            fp = futs[fut]
            parts.append(fut.result())
            print(f"  done: {os.path.basename(fp)}  ({len(parts[-1]):,} events)")
    df = pd.concat(parts, ignore_index=True)
    if cache_path:
        df.to_pickle(cache_path)
        print(f"Cached {len(df):,} events → {cache_path}")
    return df


# ── Small helpers ─────────────────────────────────────────────────────────────

def s68(x):
    """Half the 16–84 percentile span (Gaussian-equivalent sigma)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 20:
        return np.nan
    return 0.5 * (np.percentile(x, 84) - np.percentile(x, 16))


def _unit_dot(df, acols, bcols):
    A = df[list(acols)].values.astype(float)
    B = df[list(bcols)].values.astype(float)
    na = np.linalg.norm(A, axis=1)
    nb = np.linalg.norm(B, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.clip((A * B).sum(axis=1) / (na * nb), -1.0, 1.0)


def fit_median_map(d, Etrue, dmax=22.0, nbins=110, nmin=60):
    """Monotone median-regression map  edep → E_true  (the L1-optimal
    energy estimator from the deposited energy alone).
    Returns (centers, medians) for plotting + a vectorised apply function."""
    ok = np.isfinite(d) & np.isfinite(Etrue) & (d > 0)
    d, Etrue = d[ok], Etrue[ok]
    bins = np.linspace(0.0, dmax, nbins + 1)
    idx  = np.digitize(d, bins) - 1
    cs, ms = [], []
    for i in range(nbins):
        m = idx == i
        if m.sum() >= nmin:
            cs.append(0.5 * (bins[i] + bins[i + 1]))
            ms.append(np.median(Etrue[m]))
    cs = np.array(cs)
    ms = np.maximum.accumulate(np.array(ms))   # enforce monotone
    if len(cs) < 3:
        return cs, ms, lambda x: np.full_like(np.asarray(x, float), np.nan)

    def apply(x):
        x = np.asarray(x, dtype=float)
        out = np.interp(x, cs, ms, left=np.nan, right=ms[-1])
        out[~np.isfinite(x)] = np.nan
        return out
    return cs, ms, apply


def fit_linear_corr(d, Etrue, emax=22.0, nbins=44, nmin=60):
    """analyze_pairs-style: fit median upstream loss = a·E + b, invert to
    E_est = (d + b)/(1 − a).  Returns (a, b, apply)."""
    ok = np.isfinite(d) & np.isfinite(Etrue) & (d > 0)
    d, Etrue = d[ok], Etrue[ok]
    bins = np.linspace(0, emax, nbins + 1)
    idx  = np.digitize(Etrue, bins) - 1
    cs, ms = [], []
    for i in range(nbins):
        m = idx == i
        if m.sum() >= nmin and 0.5 * (bins[i] + bins[i + 1]) > 1.0:
            cs.append(0.5 * (bins[i] + bins[i + 1]))
            ms.append(np.median(Etrue[m] - d[m]))
    if len(cs) < 4:
        return None, None, lambda x: np.full_like(np.asarray(x, float), np.nan)
    a, b = np.polyfit(cs, ms, 1)

    def apply(x):
        return (np.asarray(x, dtype=float) + b) / (1.0 - a)
    return a, b, apply


def _text_page(pdf, title, elements):
    """Render a formatted text page.  elements = list of (kind, text) with
    kind ∈ {'h2', 'body', 'math', 'gap'}.  Body text is pre-wrapped."""
    fig = plt.figure(figsize=(11, 8.5))
    y = 0.965
    fig.text(0.07, y, title, fontsize=15, fontweight="bold", va="top")
    y -= 0.052
    for kind, txt in elements:
        if kind == "gap":
            y -= float(txt)
        elif kind == "h2":
            y -= 0.008
            fig.text(0.07, y, txt, fontsize=11.5, fontweight="bold",
                     va="top", color="#1a3a6b")
            y -= 0.032
        elif kind == "math":
            fig.text(0.13, y, txt, fontsize=12, va="top")
            y -= 0.042
        else:
            n = txt.count("\n") + 1
            fig.text(0.07, y, txt, fontsize=9.3, va="top", linespacing=1.38)
            y -= 0.0222 * n + 0.008
    pdf.savefig(fig)
    plt.close(fig)


def inv_mass(E1, E2, cos_th):
    return ap._inv_mass_exact(np.asarray(E1, float), np.asarray(E2, float),
                              np.asarray(cos_th, float))


def constrained_mass(E1_est, E2_est, cos_th, S=E_SUM):
    """Fix E1+E2 = S; keep only the measured sharing fraction x."""
    with np.errstate(invalid="ignore", divide="ignore"):
        x = E1_est / (E1_est + E2_est)
    return inv_mass(x * S, (1.0 - x) * S, cos_th), x


# ── Analysis ──────────────────────────────────────────────────────────────────

def run_study(df, pdf_path, rng_seed=1234):
    rng = np.random.default_rng(rng_seed)
    N   = len(df)
    et  = df["event_type"].values

    # Truth
    Et_em  = df["em_ke"].values + ME
    Et_ep  = df["ep_ke"].values + ME
    x_true = Et_em / (Et_em + Et_ep)
    cos_true  = _unit_dot(df, ("em_px", "em_py", "em_pz"),
                              ("ep_px", "ep_py", "ep_pz"))
    th_true   = np.degrees(np.arccos(cos_true))

    # Reco angles
    cos_first = _unit_dot(df, ("em_mm_px", "em_mm_py", "em_mm_pz"),
                              ("ep_mm_px", "ep_mm_py", "ep_mm_pz"))
    cos_pca   = _unit_dot(df, ("em_trk_s_px", "em_trk_s_py", "em_trk_s_pz"),
                              ("ep_trk_s_px", "ep_trk_s_py", "ep_trk_s_pz"))
    cos_nom   = _unit_dot(df, ("em_mm_gx", "em_mm_gy", "em_mm_gz"),
                              ("ep_mm_gx", "ep_mm_gy", "ep_mm_gz"))

    # Ancestry-attributed deposits (idealised attribution)
    d_em_ls  = df["em_edep_ls"].values.astype(float)
    d_ep_ls  = df["ep_edep_ls"].values.astype(float)
    d_em_all = df["em_edep_all"].values.astype(float)
    d_ep_all = df["ep_edep_all"].values.astype(float)

    # Arm-sum deposits (experimentally measurable): edep in the arm the MM
    # track entered, summed over ALL particles that deposited there.
    arm_ls  = df[[f"arm_ls_{a}"  for a in range(4)]].values.astype(float)
    arm_all = df[[f"arm_all_{a}" for a in range(4)]].values.astype(float)

    def _pick_arm(arm_table, arm_idx):
        ai  = np.asarray(arm_idx, dtype=float)
        ok  = np.isfinite(ai)
        idx = np.nan_to_num(ai, nan=0).astype(int).clip(0, 3)
        out = arm_table[np.arange(len(ai)), idx]
        out[~ok] = np.nan
        return out

    em_arm = df["em_mm_arm"].values
    ep_arm = df["ep_mm_arm"].values
    d_em_armls  = _pick_arm(arm_ls,  em_arm)
    d_ep_armls  = _pick_arm(arm_ls,  ep_arm)
    d_em_armall = _pick_arm(arm_all, em_arm)
    d_ep_armall = _pick_arm(arm_all, ep_arm)

    # ── Selection: both particles tracked in MM on different arms ───────────
    diff_arm = df["mm_diff_arm"].fillna(False).astype(bool).values
    dbl_trig = df["double_trig"].fillna(False).astype(bool).values
    sel = diff_arm & np.isfinite(cos_first)
    sel_e = sel & (d_em_armall > 0.2) & (d_ep_armall > 0.2)   # 200 keV threshold
    sig = sel_e & (et == 0)
    ipc = sel_e & (et == 1)

    print(f"\nEvents: {N:,} total | diff-arm MM {diff_arm.sum():,} "
          f"| +edep>0.2 MeV both arms {sel_e.sum():,} "
          f"(X17 {sig.sum():,}, IPC {ipc.sum():,}) "
          f"| of which double-trigger {(sel_e & dbl_trig).sum():,}")

    # ── Energy estimators (median-map fits on the selected sample) ─────────
    # Fit each estimator separately for e- and e+ (annihilation changes the
    # e+ response), on X17+IPC combined (detector property, not event type).
    fits = {}
    maps = {}
    for key, d_em, d_ep in [
        ("anc_ls",  d_em_ls,     d_ep_ls),
        ("anc_all", d_em_all,    d_ep_all),
        ("arm_ls",  d_em_armls,  d_ep_armls),
        ("arm_all", d_em_armall, d_ep_armall),
    ]:
        for p, d, Et in [("em", d_em, Et_em), ("ep", d_ep, Et_ep)]:
            cs, ms, f = fit_median_map(d[sel_e], Et[sel_e])
            fits[f"{key}_{p}"] = (cs, ms)
            maps[f"{key}_{p}"] = f
    a_em, b_em, lin_em = fit_linear_corr(d_em_ls[sel_e], Et_em[sel_e])
    a_ep, b_ep, lin_ep = fit_linear_corr(d_ep_ls[sel_e], Et_ep[sel_e])

    E_est = {
        "raw_ls":  (d_em_ls + ME,            d_ep_ls + ME),
        "lin_ls":  (lin_em(d_em_ls),         lin_ep(d_ep_ls)),
        "anc_ls":  (maps["anc_ls_em"](d_em_ls),    maps["anc_ls_ep"](d_ep_ls)),
        "anc_all": (maps["anc_all_em"](d_em_all),  maps["anc_all_ep"](d_ep_all)),
        "arm_ls":  (maps["arm_ls_em"](d_em_armls), maps["arm_ls_ep"](d_ep_armls)),
        "arm_all": (maps["arm_all_em"](d_em_armall),
                    maps["arm_all_ep"](d_ep_armall)),
    }

    # Charge-blind response map: with no magnetic field the experiment cannot
    # tell which arm holds the e+ — a single map fitted on both charges pooled
    # is what is actually implementable.
    d_pool = np.concatenate([d_em_armls[sel_e], d_ep_armls[sel_e]])
    E_pool = np.concatenate([Et_em[sel_e], Et_ep[sel_e]])
    _, _, f_cb = fit_median_map(d_pool, E_pool)
    E_est["cb_ls"] = (f_cb(d_em_armls), f_cb(d_ep_armls))

    # ── Mass ladder ─────────────────────────────────────────────────────────
    masses = {}
    masses["A truth E, truth θ"]   = inv_mass(Et_em, Et_ep, cos_true)
    masses["B truth E, θ 1st-MM"]  = inv_mass(Et_em, Et_ep, cos_first)
    masses["B2 truth E, θ PCA-smear"] = inv_mass(Et_em, Et_ep, cos_pca)
    masses["C arm-all E, truth θ"] = inv_mass(*E_est["arm_all"], cos_true)
    masses["D anc-ls E, θ 1st-MM"] = inv_mass(*E_est["anc_ls"], cos_first)
    masses["E arm-all E, θ 1st-MM"] = inv_mass(*E_est["arm_all"], cos_first)
    masses["F raw ls E, θ 1st-MM"] = inv_mass(*E_est["raw_ls"], cos_first)

    # arm-LS energies give the best sharing-fraction x (lower per-particle σ
    # than arm-all, and less compression bias)
    m_con, x_meas = constrained_mass(*E_est["arm_ls"], cos_first)
    masses["G constraint, θ 1st-MM"] = m_con
    m_con_pca, _ = constrained_mass(*E_est["arm_ls"], cos_pca)
    masses["H constraint, θ PCA-smear"] = m_con_pca
    m_cb, _ = constrained_mass(*E_est["cb_ls"], cos_first)
    masses["G2 constraint, charge-blind map"] = m_cb
    m_con_xt, _ = constrained_mass(Et_em, Et_ep, cos_first)   # perfect x
    masses["I constraint truth-x, θ 1st-MM"] = m_con_xt

    print("\n── Mass ladder (X17 signal, diff-arm + 0.2 MeV arm threshold) ──")
    print(f"{'method':38s} {'median':>8s} {'σ68':>7s} {'σ68/M':>7s}")
    ladder = {}
    for k, m in masses.items():
        v = m[sig]
        med, s = np.nanmedian(v), s68(v)
        ladder[k] = (med, s)
        print(f"{k:38s} {med:8.2f} {s:7.2f} {100*s/max(med,1e-9):6.1f}%")

    # Angle resolution numbers
    dth_first = np.degrees(np.arccos(cos_first)) - th_true
    dth_pca   = np.degrees(np.arccos(cos_pca))   - th_true
    dth_nom   = np.degrees(np.arccos(cos_nom))   - th_true
    print("\n── Opening-angle resolution (X17, σ68, deg) ──")
    for lbl, dth in [("1st-MM-hit dir", dth_first),
                     ("PCA fit (0.5 mm smear)", dth_pca),
                     ("target-centre line", dth_nom)]:
        print(f"  {lbl:28s} {s68(dth[sig]):6.2f}")

    # Per-particle energy resolution after correction
    print("\n── Per-particle energy resolution σ68(E_est−E_true) [MeV] "
          "(X17+IPC, sel) ──")
    for key in ["raw_ls", "lin_ls", "anc_ls", "anc_all", "arm_ls", "arm_all"]:
        r_em = E_est[key][0][sel_e] - Et_em[sel_e]
        r_ep = E_est[key][1][sel_e] - Et_ep[sel_e]
        print(f"  {key:8s}  e-: {s68(r_em):5.2f}   e+: {s68(r_ep):5.2f}")

    # ── Photostatistics scenarios ────────────────────────────────────────────
    # Smear the arm-sum deposits with σ_E = k·√E (E in MeV) — stochastic
    # photoelectron term; k=0.10 ⇒ 3.2% at 10 MeV (≈100 p.e./MeV).
    ks = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50]
    photo_direct, photo_con = [], []
    for k in ks:
        if k == 0:
            dem, dep = d_em_armls, d_ep_armls
        else:
            dem = d_em_armls + rng.normal(0, 1, N) * k * np.sqrt(
                np.clip(d_em_armls, 0, None))
            dep = d_ep_armls + rng.normal(0, 1, N) * k * np.sqrt(
                np.clip(d_ep_armls, 0, None))
        _, _, f_em = fit_median_map(dem[sel_e], Et_em[sel_e])
        _, _, f_ep = fit_median_map(dep[sel_e], Et_ep[sel_e])
        E1, E2 = f_em(dem), f_ep(dep)
        m_dir = inv_mass(E1, E2, cos_first)
        m_c, _ = constrained_mass(E1, E2, cos_first)
        photo_direct.append((np.nanmedian(m_dir[sig]), s68(m_dir[sig])))
        photo_con.append((np.nanmedian(m_c[sig]), s68(m_c[sig])))
    print("\n── Photostatistics scan (arm-ls estimator, X17 σ68 [MeV]) ──")
    print(f"{'k [√MeV]':>9s} {'σ/E @10MeV':>11s} {'direct':>8s} {'constraint':>11s}")
    for k, (md, sd), (mc, sc) in zip(ks, photo_direct, photo_con):
        print(f"{k:9.2f} {100*k/np.sqrt(10):10.1f}% {sd:8.2f} {sc:11.2f}")

    # ── Analytic budget ──────────────────────────────────────────────────────
    sE_em = s68((E_est["arm_all"][0] - Et_em)[sig]) / np.nanmedian(Et_em[sig])
    sE_ep = s68((E_est["arm_all"][1] - Et_ep)[sig]) / np.nanmedian(Et_ep[sig])
    sth   = np.radians(s68(dth_first[sig]))
    th0   = np.radians(np.nanmedian(th_true[sig]))
    ang_term = 0.5 / np.tan(th0 / 2) * sth
    en_term  = 0.5 * np.hypot(sE_em, sE_ep)
    print("\n── Analytic budget  σM/M ≈ ½σE1/E1 ⊕ ½σE2/E2 ⊕ ½cot(θ/2)σθ ──")
    print(f"  σE/E (e-): {100*sE_em:.1f}%  (e+): {100*sE_ep:.1f}%  "
          f"σθ: {np.degrees(sth):.2f}°  θ̄: {np.degrees(th0):.1f}°")
    print(f"  energy term: {100*en_term:.1f}%   angular term: {100*ang_term:.1f}%"
          f"   total: {100*np.hypot(en_term, ang_term):.1f}%"
          f"   (× {M_X17:.1f} MeV = {M_X17*np.hypot(en_term, ang_term):.2f} MeV)")

    # ──────────────────────────────────────────────────────────────────────────
    # Plots
    # ──────────────────────────────────────────────────────────────────────────
    mass_bins = np.linspace(0, 25, 126)
    mc = 0.5 * (mass_bins[:-1] + mass_bins[1:])

    # Numbers quoted in the explanation pages
    _s_raw  = s68((E_est["raw_ls"][0] - Et_em)[sel_e])
    _s_lin  = s68((E_est["lin_ls"][0] - Et_em)[sel_e])
    _s_map  = s68((E_est["arm_ls"][0] - Et_em)[sel_e])
    _sB  = ladder["B truth E, θ 1st-MM"][1]
    _sE  = ladder["E arm-all E, θ 1st-MM"][1]
    _sG, _mG  = ladder["G constraint, θ 1st-MM"][1], ladder["G constraint, θ 1st-MM"][0]
    _sG2 = ladder["G2 constraint, charge-blind map"][1]
    _sF  = ladder["F raw ls E, θ 1st-MM"][1]

    with PdfPages(pdf_path) as pdf:

        # ── Explanation page 1: what this study does ────────────────────────
        _text_page(pdf, "1.  What this study does", [
            ("body",
             "Goal: determine how well the e⁺e⁻ invariant mass can be reconstructed with the scintillator\n"
             "calorimetry, identify what limits the resolution, and define the algorithm the experiment should\n"
             "implement.  Input: 400k Geant4 events (4 local job files, old geometry) — correlated pairs from the\n"
             "20.58 MeV ⁴He* transition, 50% X17 (m = 16.8 MeV) / 50% IPC (dN/dM ∝ 1/M), uniform vertex in the\n"
             "He-3 gas.  Selection: both tracks seen in the Micromegas on two different arms, ≥0.2 MeV deposited\n"
             f"in both arms  (N = {sig.sum():,} X17, {ipc.sum():,} IPC)."),
            ("h2", "The three measured ingredients"),
            ("math",
             r"$M^2 \;=\; 2m_e^2 \;+\; 2\,\left(E_1 E_2 \,-\, p_1 p_2 \cos\theta\right)$"),
            ("body",
             f"•  Opening angle θ from the two MM tracks.  Resolution σ68 ≈ {s68(dth_first[sig]):.1f}°, dominated by multiple\n"
             "    scattering in the target wall and entrance materials — NOT by MM spatial resolution (the PCA track\n"
             "    fit and the ideal first-hit direction give nearly identical residuals, page 6).\n"
             "•  Lepton energies from calorimetry.  Two attributions are compared: summing deposits per ARM (the\n"
             "    only thing the experiment can measure) vs per Geant4 ancestry (idealised).  They agree closely —\n"
             "    cross-arm shower leakage is negligible — so per-arm summing is a valid measurement.\n"
             "•  An energy response correction, derived from simulation (page 4 / section 3)."),
            ("h2", "Why the calorimeter alone is weak — and what rescues it"),
            ("body",
             "The LS stack is ≈ 0.1 radiation lengths: it range-stops electrons but most bremsstrahlung escapes,\n"
             f"so the per-particle energy resolution is only σ_E/E ≈ 20–30% ({_s_map:.1f} MeV σ68 for e⁻ after correction).\n"
             f"A direct two-energy mass reconstruction therefore gives σ68 ≈ {_sE:.1f} MeV.  The rescue is kinematic:\n"
             "the total pair energy is fixed by the transition (page 2), which reduces the calorimeter's job to\n"
             f"measuring an energy-sharing RATIO.  This brings the resolution to σ68 ≈ {_sG:.1f} MeV, close to the\n"
             f"angular floor of {_sB:.2f} MeV.  Photostatistics hardly matters (page 10): leakage fluctuations dominate,\n"
             "so light-collection design should be driven by thresholds and timing, not energy resolution."),
            ("h2", "Reading the rest of this document"),
            ("body",
             "Pages 4–5: calorimeter response maps and per-particle energy estimators;  page 6: angle residuals.\n"
             "Page 7: the error-budget 'ladder' — truth quantities are swapped for reconstructed ones one at a\n"
             "time, so each error source's contribution to the mass width is isolated.  Page 8: the energy-sharing\n"
             "fraction.  Page 9: X17 vs IPC separation and mass-window ROC.  Page 11: summary numbers."),
        ])

        # ── Explanation page 2: the constraint ──────────────────────────────
        _text_page(pdf, "2.  The 20.58 MeV constraint — use and bias", [
            ("h2", "Where it comes from"),
            ("body",
             "Slow-neutron capture, ³He + n → ⁴He*, produces the 20.58 MeV state essentially at rest (the n_TOF\n"
             "neutron contributes ≲ keV in the CM).  The transition energy is shared between the e⁺e⁻ pair and the\n"
             "⁴He recoil — but the ⁴He is 3727 MeV heavy, so the recoil T = p²/2M(⁴He) is tiny:  18.8 keV for X17\n"
             "(p = 11.9 MeV), and for IPC it varies with the virtual-photon mass from 56 keV (Mee → 2mₑ, p ≈ 20.6)\n"
             "down to 0 (Mee at the endpoint).  The pair energy sum is therefore NOT exactly equal for the two\n"
             "modes — but the full spread is < 0.06 MeV, and the X17-vs-IPC difference (≤ 38 keV) propagates to a\n"
             "constrained-mass shift of only ~30 keV: a factor ~50 below the resolution.  (This simulation does not\n"
             "model the recoil at all — S = 20.58 MeV exactly for both modes.)  Hence, event by event:"),
            ("math",
             r"$S \;\equiv\; E(e^+) + E(e^-) \;=\; 20.58\ \mathrm{MeV} \;-\; T_{rec}(M_{ee})\,,\qquad T_{rec} < 57\ \mathrm{keV}$"),
            ("h2", "How it is used"),
            ("body",
             "Write E₁ = xS and E₂ = (1−x)S.  The mass becomes a function of just two variables, the sharing\n"
             "fraction x and the opening angle θ (ultra-relativistic form shown for intuition; the exact formula\n"
             "with mₑ is used in the code):"),
            ("math",
             r"$M(x,\theta) \;\approx\; S\,\sin(\theta/2)\,\sqrt{4x(1-x)}\,,\qquad \hat{x} \;=\; \hat{E}_1/(\hat{E}_1+\hat{E}_2)$"),
            ("body",
             "Three things make this powerful:\n"
             "•  ∂M/∂x ∝ (1−2x) = 0 at x = ½ — first-order INSENSITIVITY to calorimeter error, exactly where the\n"
             "    X17 population is densest (x is uniform on [0.21, 0.79] for the β = 0.578 X17 boost).\n"
             "•  x̂ is a ratio: correlated calibration errors (common gain drifts, shared response biases) cancel.\n"
             "•  The absolute pair energy scale is set by S, which is known — the calorimeter's poor absolute\n"
             f"    resolution no longer enters.  Measured effect: σ68 {_sE:.2f} → {_sG:.2f} MeV  (truth-x limit: {ladder['I constraint truth-x, θ 1st-MM'][1]:.2f})."),
            ("h2", "Does it introduce bias?  Three separate questions"),
            ("body",
             "1.  The constraint itself: NO bias, provided the event really is a pair from this transition.  S is\n"
             "     then the true sum, and unbiased (x̂, θ̂) give an unbiased mass.  Crucially the IPC background\n"
             "     comes from the SAME transition and satisfies the same constraint — the signal/background\n"
             "     comparison stays fair, and the IPC spectrum remains broad (page 9).\n"
             f"2.  Estimator bias: the median-map energy correction compresses x̂ toward ½ (regression to the\n"
             f"     mean), and M is maximal at x = ½ — so the X17 peak reconstructs ≈ {_mG-M_X17:+.1f} MeV high ({_mG:.1f} vs 16.8).\n"
             "     This is a smooth, simulation-known SCALE shift, not a distortion of the ordering: fold it into\n"
             "     the fit template or calibrate the scale on simulation.  Separation power is unaffected.\n"
             "3.  Events that violate the constraint (random coincidences, pile-up, cosmics): for them the\n"
             "     constrained mass is a deformed variable with no particular peak.  Better: the measured raw sum\n"
             "     Ê₁+Ê₂ is an independent observable — requiring it to be compatible with 20.58 MeV rejects\n"
             "     non-transition backgrounds BEFORE the constraint is applied, and its data distribution is an\n"
             "     in-situ standard candle for validating the simulation-derived energy correction.\n"
             "\n"
             "A free bonus: constrained M ≤ S by construction, so the unphysical tail above the kinematic\n"
             "endpoint that the direct method produces (page 7, red) disappears automatically."),
        ])

        # ── Explanation page 3: implementation recipe ───────────────────────
        _text_page(pdf, "3.  Implementing it in the experiment", [
            ("body",
             "Assume per-cell energy calibrations are known (source/cosmic calibration of each LS cell, trigger\n"
             "scintillator and back scintillator, in MeV electron-equivalent).  The pipeline per event is:"),
            ("h2", "Event-level algorithm"),
            ("body",
             "1.  Select pairs: MM track on two different arms, trigger coincidence.\n"
             "2.  Arm energy: dᵢ = calibrated sum of the LS cells of arm i.  (LS-only gives the best sharing\n"
             "     fraction; trigger-scint and back-scint sums are kept for the Ê₁+Ê₂ consistency check.)\n"
             "3.  Response correction: Êᵢ = f(dᵢ), where f is the monotone median map  f(d) = median(E_true | d)\n"
             "     DERIVED FROM GEANT4.  One charge-blind map (no B-field → arm charge unknown):  costs\n"
             f"     essentially nothing, σ68 = {_sG2:.2f} vs {_sG:.2f} MeV with per-charge truth maps.\n"
             "4.  Opening angle θ from the two MM track directions.\n"
             "5.  x̂ = Ê₁/(Ê₁+Ê₂);   M̂ = exact mass formula with E₁ = x̂S, E₂ = (1−x̂)S, S = 20.58 MeV.\n"
             "6.  Inference: fit the data M̂ spectrum with X17(m) + IPC templates produced by the SAME simulation\n"
             "     chain — the estimator scale shift of point 2 (page 2) lives inside the template, so it does not\n"
             "     bias the extracted mass or yield."),
            ("h2", "Why a median map and not the linear loss fit"),
            ("body",
             "The linear inversion E = (d + b)/(1−a) removes the average upstream loss but AMPLIFIES every\n"
             f"fluctuation by 1/(1−a):  per-particle σ68 = {_s_lin:.1f} MeV, vs {_s_raw:.1f} raw and {_s_map:.1f} for the median map.\n"
             "The median map is the L1-optimal point estimate of E given d, and being monotone it is trivially\n"
             "invertible and stable.  Its one side effect is the regression-to-the-mean compression discussed on\n"
             "page 2 — handled by template fitting."),
            ("h2", "What must come from simulation, and how to validate it in-situ"),
            ("body",
             "Yes — the loss before the calorimeter (gas window, drift gas, mesh, PCB stack, air gaps, wrapping,\n"
             "trigger scintillator, CFRP/Al liners) must be calibrated via simulation: it cannot be measured event\n"
             "by event.  The map f bundles upstream loss, longitudinal/lateral leakage and e⁺ annihilation into\n"
             "one response function.  Practical consequences:\n"
             "•  Fold the measured light response (Birks quenching, collection non-uniformity, resolution) into\n"
             "    the simulation BEFORE deriving f — this study uses bare edep; the photostatistics scan (page 10)\n"
             "    shows the readout term is subdominant, but quenching shifts the scale and must be in the model.\n"
             "•  Re-derive f and the templates after EVERY geometry change (same script, --rebuild).\n"
             "•  In-situ validation handles, in order of power:\n"
             "      a.  the measured-sum spectrum Ê₁+Ê₂ must center on 20.58 MeV — a shift directly measures a\n"
             "           response-correction error, per arm pair, using the data themselves;\n"
             "      b.  the IPC-dominated region (low constrained M) is a known QED spectrum — a shape standard\n"
             "           candle for the full reconstruction chain;\n"
             "      c.  single-particle calibration runs (known energy/angle, --single mode) per arm."),
        ])

        # ── Page 1: energy response maps ────────────────────────────────────
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle("Calorimeter response: deposited energy vs true E_total "
                     "(X17+IPC, diff-arm selection)", fontsize=12)
        for axc, (ttl, d, Et, key, p) in zip(axes.flat, [
            ("e⁻ ancestry LS edep",      d_em_ls,     Et_em, "anc_ls",  "em"),
            ("e⁺ ancestry LS edep",      d_ep_ls,     Et_ep, "anc_ls",  "ep"),
            ("e⁻ arm-sum all-scint edep", d_em_armall, Et_em, "arm_all", "em"),
            ("e⁺ arm-sum all-scint edep", d_ep_armall, Et_ep, "arm_all", "ep"),
        ]):
            ok = sel_e & np.isfinite(d) & (d > 0)
            h, xe, ye = np.histogram2d(Et[ok], d[ok],
                                       bins=[np.linspace(0, 22, 88)] * 2)
            rs = h.sum(axis=1, keepdims=True).clip(1)
            axc.pcolormesh(xe, ye, (h / rs).T, cmap="viridis")
            axc.plot([0, 22], [0, 22], "w--", lw=1)
            cs, ms = fits[f"{key}_{p}"]
            if len(cs):
                axc.plot(ms, cs, "r-", lw=1.5, label="median map (inverted)")
            axc.set_xlabel("True E_total [MeV]")
            axc.set_ylabel("edep [MeV]")
            axc.set_title(ttl)
            axc.legend(fontsize=8, loc="upper left")
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # ── Page 2: per-particle energy resolution vs E_true ────────────────
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
        ebins = np.linspace(1, 21, 21)
        ecen  = 0.5 * (ebins[:-1] + ebins[1:])
        for axc, p, Et in [(axes[0], 0, Et_em), (axes[1], 1, Et_ep)]:
            for key, colr in [("raw_ls", "#999999"), ("lin_ls", "#d62728"),
                              ("anc_ls", "#1f77b4"), ("arm_ls", "#2ca02c"),
                              ("arm_all", "#9467bd")]:
                res = E_est[key][p] - Et
                ss = []
                for i in range(len(ecen)):
                    m = sel_e & (Et >= ebins[i]) & (Et < ebins[i + 1])
                    ss.append(s68(res[m]))
                axc.plot(ecen, ss, "-o", ms=3, color=colr, label=key)
            axc.set_xlabel("True E_total [MeV]")
            axc.set_ylabel("σ68(E_est − E_true) [MeV]")
            axc.set_title("e⁻" if p == 0 else "e⁺")
            axc.grid(alpha=0.3); axc.legend(fontsize=8)
        fig.suptitle("Per-particle energy resolution by estimator", fontsize=12)
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # ── Page 3: angle resolution recap ──────────────────────────────────
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
        ax = axes[0]
        bins = np.linspace(-25, 25, 151)
        for dth, lbl, colr in [(dth_first, "1st-MM-hit dir", "#000000"),
                               (dth_pca, "PCA fit (0.5 mm smear)", "#2ca02c"),
                               (dth_nom, "target centre → 1st hit", "#9467bd")]:
            v = dth[sig]
            ax.hist(v[np.isfinite(v)], bins=bins, histtype="step", lw=1.8,
                    color=colr, density=True,
                    label=f"{lbl}  σ68={s68(v):.2f}°")
        ax.set_xlabel("θ_reco − θ_true [deg]"); ax.set_ylabel("density")
        ax.set_title("Opening-angle residuals (X17, diff-arm)")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

        ax = axes[1]
        kemin = np.minimum(df["em_ke"].values, df["ep_ke"].values)
        kb = np.linspace(0, 11, 12)
        kc = 0.5 * (kb[:-1] + kb[1:])
        for dth, lbl, colr in [(dth_first, "1st-MM-hit dir", "#000000"),
                               (dth_pca, "PCA fit", "#2ca02c")]:
            ss = [s68(dth[sig & (kemin >= kb[i]) & (kemin < kb[i + 1])])
                  for i in range(len(kc))]
            ax.plot(kc, ss, "-o", ms=4, color=colr, label=lbl)
        ax.set_xlabel("min(KE e⁻, KE e⁺) [MeV]")
        ax.set_ylabel("σ68(Δθ) [deg]")
        ax.set_title("Angle resolution vs softer-particle KE")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # ── Page 4: mass ladder ─────────────────────────────────────────────
        fig, axes = plt.subplots(1, 2, figsize=(15, 5.8))
        order = ["B truth E, θ 1st-MM", "C arm-all E, truth θ",
                 "E arm-all E, θ 1st-MM", "F raw ls E, θ 1st-MM"]
        colors = ["#2ca02c", "#1f77b4", "#d62728", "#999999"]
        ax = axes[0]
        for k, colr in zip(order, colors):
            v = masses[k][sig]
            ax.hist(v[np.isfinite(v)], bins=mass_bins, histtype="step", lw=1.8,
                    color=colr, density=True,
                    label=f"{k}  σ68={ladder[k][1]:.2f}")
        ax.axvline(M_X17, color="red", ls=":", lw=1.2)
        ax.set_xlabel("M(e⁺e⁻) [MeV]"); ax.set_ylabel("density")
        ax.set_title("Direct reconstruction — error-budget ladder (X17)")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

        ax = axes[1]
        for k, colr in [("E arm-all E, θ 1st-MM", "#d62728"),
                        ("G constraint, θ 1st-MM", "#1f77b4"),
                        ("G2 constraint, charge-blind map", "#ff7f0e"),
                        ("H constraint, θ PCA-smear", "#17becf"),
                        ("I constraint truth-x, θ 1st-MM", "#2ca02c")]:
            v = masses[k][sig]
            ax.hist(v[np.isfinite(v)], bins=mass_bins, histtype="step", lw=1.8,
                    color=colr, density=True,
                    label=f"{k}  σ68={ladder[k][1]:.2f}")
        ax.axvline(M_X17, color="red", ls=":", lw=1.2)
        ax.set_xlabel("M(e⁺e⁻) [MeV]"); ax.set_ylabel("density")
        ax.set_title("E₁+E₂ = 20.58 MeV constraint (X17)")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        fig.suptitle("Invariant-mass reconstruction ladder", fontsize=12)
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # ── Page 5: energy sharing x ────────────────────────────────────────
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
        ax = axes[0]
        xb = np.linspace(0, 1, 81)
        ax.hist(x_true[sig], bins=xb, histtype="step", color="#2ca02c", lw=1.8,
                density=True, label="x_true (X17)")
        ax.hist(x_meas[sig & np.isfinite(x_meas)], bins=xb, histtype="step",
                color="#d62728", lw=1.8, density=True,
                label="x_meas (arm-ls; spikes = degenerate\nlow-edep ends of the median map)")
        ax.hist(x_true[ipc], bins=xb, histtype="step", color="#1f77b4", lw=1.4,
                ls="--", density=True, label="x_true (IPC)")
        ax.set_xlabel("x = E(e⁻)/(E(e⁻)+E(e⁺))"); ax.set_ylabel("density")
        ax.set_title("Energy-sharing fraction"); ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

        ax = axes[1]
        dx = (x_meas - x_true)[sig]
        ax.hist(dx[np.isfinite(dx)], bins=np.linspace(-0.4, 0.4, 121),
                histtype="step", color="#d62728", lw=1.8, density=True,
                label=f"σ68 = {s68(dx):.3f}")
        ax.set_xlabel("x_meas − x_true"); ax.set_ylabel("density")
        ax.set_title("Sharing-fraction resolution (X17)\n"
                     "M is 1st-order insensitive to x near x = ½")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # ── Page 6: X17 vs IPC separation ───────────────────────────────────
        fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
        roc_methods = [("E arm-all E, θ 1st-MM", "direct, arm-all", "#d62728"),
                       ("G constraint, θ 1st-MM", "constraint", "#1f77b4"),
                       ("H constraint, θ PCA-smear", "constraint+PCA", "#17becf")]
        for axc, (k, lbl, colr) in zip(axes[:2], roc_methods[:2]):
            for mask, elbl, ecol, els in [(sig, "X17", "#e84040", "-"),
                                          (ipc, "IPC", "#4a90d9", "--")]:
                v = masses[k][mask]
                axc.hist(v[np.isfinite(v)], bins=mass_bins, histtype="step",
                         lw=1.8, ls=els, color=ecol, density=True, label=elbl)
            axc.axvline(M_X17, color="red", ls=":", lw=1.2)
            axc.set_xlabel("M(e⁺e⁻) [MeV]"); axc.set_ylabel("density")
            axc.set_title(f"X17 vs IPC — {lbl}")
            axc.legend(fontsize=9); axc.grid(alpha=0.3)

        ax = axes[2]
        for k, lbl, colr in roc_methods:
            ms_, mb_ = masses[k][sig], masses[k][ipc]
            ms_, mb_ = ms_[np.isfinite(ms_)], mb_[np.isfinite(mb_)]
            med = np.median(ms_)
            ws = np.linspace(0.05, 12, 200)
            eff_s = [(np.abs(ms_ - med) < w).mean() for w in ws]
            eff_b = [(np.abs(mb_ - med) < w).mean() for w in ws]
            ax.plot(eff_s, eff_b, color=colr, lw=2, label=lbl)
            i90 = np.searchsorted(eff_s, 0.90)
            if i90 < len(ws):
                print(f"  ROC {lbl:18s}: IPC eff at 90% X17 eff = "
                      f"{eff_b[i90]*100:.1f}%  (|ΔM| < {ws[i90]:.2f} MeV)")
        ax.set_xlabel("X17 signal efficiency (mass window)")
        ax.set_ylabel("IPC efficiency (contamination)")
        ax.set_yscale("log")
        ax.set_title("Mass-window ROC\n(generator IPC: dN/dM ∝ 1/M)")
        ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # ── Page 7: photostatistics scan ────────────────────────────────────
        fig, ax = plt.subplots(figsize=(9, 5.5))
        ax.plot(ks, [s for _, s in photo_direct], "-o", color="#d62728",
                label="direct (arm-ls E + θ 1st-MM)")
        ax.plot(ks, [s for _, s in photo_con], "-o", color="#1f77b4",
                label="constrained (E₁+E₂ = 20.58 MeV)")
        ax2 = ax.secondary_xaxis(
            "top", functions=(lambda k: 100 * k / np.sqrt(10),
                              lambda p: p * np.sqrt(10) / 100))
        ax2.set_xlabel("σ_E/E at 10 MeV [%]")
        ax.set_xlabel("photostatistics k  (σ_E = k·√E, E in MeV)")
        ax.set_ylabel("X17 mass σ68 [MeV]")
        ax.set_title("Impact of light-collection statistics on mass resolution")
        ax.legend(fontsize=9); ax.grid(alpha=0.3)
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # ── Page 8: summary text ────────────────────────────────────────────
        fig = plt.figure(figsize=(11, 8.5))
        lines = [
            "INVARIANT-MASS RECONSTRUCTION SUMMARY",
            "",
            f"Selection: both tracks in MM on different arms, >0.2 MeV in both "
            f"arms  (N_X17 = {sig.sum():,}, N_IPC = {ipc.sum():,})",
            "",
            "Mass ladder (X17): median ± σ68 [MeV]",
        ] + [f"   {k:38s} {ladder[k][0]:6.2f} ± {ladder[k][1]:5.2f}"
             for k in ladder] + [
            "",
            f"Angle σ68: 1st-MM {s68(dth_first[sig]):.2f}° | "
            f"PCA-smear {s68(dth_pca[sig]):.2f}° | "
            f"target-line {s68(dth_nom[sig]):.2f}°",
            f"Energy σ68/E (arm-all, corrected): e⁻ {100*sE_em:.1f}%, "
            f"e⁺ {100*sE_ep:.1f}%",
            "",
            "Analytic: σM/M ≈ ½σE1/E1 ⊕ ½σE2/E2 ⊕ ½cot(θ/2)·σθ",
            f"   energy term {100*en_term:.1f}% | angle term {100*ang_term:.1f}% "
            f"| total {100*np.hypot(en_term, ang_term):.1f}% "
            f"→ {M_X17*np.hypot(en_term, ang_term):.2f} MeV",
            "",
            "Implementation recipe:",
            " 1. Energy per particle = arm-summed LS edep in the arm the MM",
            "    track points to (LS-only beats adding PlScint+BS for x).",
            " 2. Correct with the monotone median map edep → E_total (per charge).",
            " 3. Opening angle from MM track directions.",
            " 4. Apply E₁+E₂ = 20.58 MeV: M = M(xS, (1−x)S, θ) with x = Ê₁/(Ê₁+Ê₂).",
            " 5. Calibrate the mass scale on simulation (median-map compression",
            "    biases the peak high by ~0.5–1 MeV; width is what matters).",
            "",
            "Caveats: no optical transport/Birks (photostat scan approximates it);",
            "ancestry vs arm-sum quantifies attribution cross-talk; old geometry;",
            "1st-MM-hit direction is slightly optimistic (use PCA row as realistic).",
        ]
        fig.text(0.06, 0.96, "\n".join(lines), va="top", family="monospace",
                 fontsize=9)
        pdf.savefig(fig); plt.close(fig)

    print(f"\nPDF written: {pdf_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    pa = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    pa.add_argument("inputs", nargs="+", help="ROOT files (or dir)")
    pa.add_argument("-o", "--output", default="invmass_study.pdf")
    pa.add_argument("--cache", default=None,
                    help="per-event table pickle (default: alongside output)")
    pa.add_argument("--rebuild", action="store_true",
                    help="ignore existing cache")
    pa.add_argument("--cache-only", action="store_true",
                    help="build the cache and exit")
    pa.add_argument("--workers", type=int, default=4)
    args = pa.parse_args()

    files = ap.collect_files(args.inputs, "*.root")
    if not files:
        sys.exit("No input files found")
    cache = args.cache or os.path.splitext(args.output)[0] + "_perevent.pkl"

    df = build_table(files, cache, rebuild=args.rebuild, workers=args.workers)
    if args.cache_only:
        return
    run_study(df, args.output)


if __name__ == "__main__":
    main()
