#!/usr/bin/env python3
"""
analyze_pairs.py — MX17 Geant4 X17+IPC pair simulation analysis
=================================================================
Reads HitTree + EventTree from mx17_full_sim output and answers:

  Acceptance    — do e± reach each detector layer?
  Calorimetry   — does the LS stack reconstruct a useful invariant mass?
  Calo QA       — LS energy containment, stopping fraction
  Asymmetry     — are low-energy particles in asymmetric pairs stopped?
  Pointing      — do MM tracks point back to the He-3 vertex?
  Angle reco    — how well can we reconstruct the true opening angle?

Processes multiple unmerged ROOT files without loading everything into memory.
HitTree is read in 300k-row chunks; only the columns needed are loaded.

Usage:
    python analyze_pairs.py /eos/.../pairs/x17_ipc_pairs_job*.root
    python analyze_pairs.py /eos/.../pairs/ -o analysis.pdf --max-files 20
"""

import argparse
import gc
import glob
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import uproot
import pandas as pd
from tqdm import tqdm

# ── Constants ─────────────────────────────────────────────────────────────────
ME_MEV         = 0.511
# A primary particle whose minimum recorded KE (across all scored-volume steps) is
# below this threshold is considered to have stopped in that volume.
STOP_KE_THRESH = 0.1  # MeV
TRIGGER_LAYERS = frozenset(["PlasticScint", "LiqScint_1"])
LS_LAYERS      = frozenset(["LiqScint_1", "LiqScint_2"])
BS_LAYERS      = frozenset(["BackScintL", "BackScintR"])

def _inv_mass_exact(E1, E2, cos_th):
    """Exact two-particle invariant mass including electron/positron rest mass.

    M² = 2·mₑ² + 2·(E1·E2 − |p1|·|p2|·cosθ)    where |p_i| = sqrt(E_i² − mₑ²)

    Using the ultra-relativistic limit |p| ≈ E underestimates M at low energies
    (5–10% for E ≈ 1.5 MeV, <0.2% for E > 5 MeV).
    """
    p1 = np.sqrt(np.maximum(E1**2 - ME_MEV**2, 0.0))
    p2 = np.sqrt(np.maximum(E2**2 - ME_MEV**2, 0.0))
    M2 = 2*ME_MEV**2 + 2*(E1*E2 - p1*p2*cos_th)
    return np.sqrt(np.maximum(M2, 0.0))


ETYPE_LABEL = {0: "X17 signal", 1: "IPC background"}
ETYPE_COLOR = {0: "#e84040",    1: "#4a90d9"}
ETYPE_LS    = {0: "-",          1: "--"}

# ── Angle-reconstruction study constants ─────────────────────────────────────
# Proxy for MM cluster position resolution applied to drift hit positions
# before the smeared PCA fit (per-step smearing, crude digitization stand-in).
MM_HIT_SMEAR_MM = 0.5
# Trigger-scintillator position resolution: ~10 cm segmentation → σ ≈ 10/√12 cm
SCINT_SMEAR_MM  = 30.0

# Per-particle direction estimators: ψ = space angle between two directions.
# All are compared to the TRUTH direction at the vertex unless noted.
PSI_ESTIMATORS = ["first", "last", "drift", "fit", "fit_vs_first",
                  "vline", "vcen", "nomline", "scint"]
PSI_LABEL = {
    "first":        "dir @ 1st MM hit (best measurable)",
    "last":         "dir @ last MM hit",
    "drift":        "1st vs last MM dir (drift-gas MS)",
    "fit":          "PCA fit of drift hits",
    "fit_vs_first": "PCA fit vs 1st-hit dir (fit error only)",
    "vline":        "true vtx → 1st MM hit",
    "vcen":         "true vtx → drift edep centroid",
    "nomline":      "target centre → 1st MM hit",
    "scint":        f"true vtx → trig-scint hit (σ={SCINT_SMEAR_MM/10:.0f} cm)",
}
PSI_COLOR = {
    "first":        "#000000",
    "last":         "#7f7f7f",
    "drift":        "#8c564b",
    "fit":          "#2ca02c",
    "fit_vs_first": "#98df8a",
    "vline":        "#1f77b4",
    "vcen":         "#17becf",
    "nomline":      "#9467bd",
    "scint":        "#ff7f0e",
}

# Opening-angle reconstruction methods (need both particles)
OPEN_METHODS = ["first", "fit", "fits", "vline", "vcen", "nomline", "scint"]
METH_LABEL = {
    "first":   "current: dir @ 1st MM hit",
    "fit":     "PCA drift-track fit",
    "fits":    f"PCA fit, {MM_HIT_SMEAR_MM:.1f} mm hit smear",
    "vline":   "true vtx → 1st MM hit",
    "vcen":    "true vtx → drift centroid",
    "nomline": "target centre → 1st MM hit",
    "scint":   f"true vtx → scint hit (σ={SCINT_SMEAR_MM/10:.0f} cm)",
}
METH_COLOR = {
    "first":   "#000000",
    "fit":     "#2ca02c",
    "fits":    "#98df8a",
    "vline":   "#1f77b4",
    "vcen":    "#17becf",
    "nomline": "#9467bd",
    "scint":   "#ff7f0e",
}


def _angle_deg(A, B):
    """Space angle [deg] between row vectors of A and B (need not be unit).
    NaN rows propagate to NaN."""
    na = np.linalg.norm(A, axis=1)
    nb = np.linalg.norm(B, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (A * B).sum(axis=1) / (na * nb)
        return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def _arr(m, cols):
    """Extract columns as float array; all-NaN if any column is missing."""
    if all(c in m.columns for c in cols):
        return m[list(cols)].values.astype(float)
    return np.full((len(m), len(cols)), np.nan)

HIT_COLS = ["eventID", "trackID", "parentID", "armID",
            "detType", "particle", "edep", "ke",
            "gx", "gy", "gz", "px", "py", "pz"]

SCORED_LAYERS = ["DriftGas", "PlasticScint", "LiqScint_1",
                 "LiqScint_2", "BackScintL", "BackScintR"]
LAYER_LABELS  = {
    "DriftGas":    "MM DriftGas",
    "PlasticScint":"Trigger scint.",
    "LiqScint_1":  "Liq. scint. 1",
    "LiqScint_2":  "Liq. scint. 2",
    "BackScintL":  "Back scint. L",
    "BackScintR":  "Back scint. R",
}
LAYER_COLOR = {
    "DriftGas":    "#2ca02c",
    "PlasticScint":"#ff7f0e",
    "LiqScint_1":  "#6a1d8a",
    "LiqScint_2":  "#9b59b6",
    "BackScintL":  "#d62728",
    "BackScintR":  "#e57373",
}


# ── String decoding ───────────────────────────────────────────────────────────

def _decode_col(series):
    """Convert uproot char-array column (bytes or str) to plain str."""
    if series.empty:
        return series
    sample = series.iloc[0]
    if isinstance(sample, bytes):
        return series.apply(lambda x: x.rstrip(b"\x00").decode("ascii", errors="ignore"))
    return series.str.rstrip("\x00") if hasattr(series, "str") else series.astype(str)


# ── Hit processing (per chunk) ────────────────────────────────────────────────

def _drift_track_quantities(mm, rng):
    """Per-event drift-track quantities from one primary's DriftGas hits
    (restricted to a single arm by the caller).

    Returns a DataFrame indexed by eventID with columns:
      n               number of drift hits
      l_px,l_py,l_pz  momentum direction at the LAST drift hit
      c_gx,c_gy,c_gz  edep-weighted centroid of the drift hit positions
      f_px,f_py,f_pz  PCA straight-line fit direction (NaN if n < 3)
      s_px,s_py,s_pz  PCA fit after smearing positions by MM_HIT_SMEAR_MM
    PCA directions are oriented along the first-hit momentum direction.
    The covariance/eigenvector computation is batched over all events with
    np.bincount + np.linalg.eigh, so cost is O(hits), not O(events·hits).
    """
    g = mm.groupby("eventID")
    out = pd.DataFrame(index=g.size().index)
    out["n"] = g.size().astype(float)
    out[["l_px", "l_py", "l_pz"]] = g.last()[["px", "py", "pz"]]

    eids, inv = np.unique(mm["eventID"].values, return_inverse=True)
    X = mm[["gx", "gy", "gz"]].values.astype(float)
    n = np.bincount(inv).astype(float)

    # edep-weighted centroid
    w = mm["edep"].values.astype(float)
    wsum = np.bincount(inv, w).clip(1e-30)
    cen = np.stack([np.bincount(inv, w * X[:, k]) for k in range(3)], axis=1)
    out[["c_gx", "c_gy", "c_gz"]] = cen / wsum[:, None]

    d_first = g.first()[["px", "py", "pz"]].values.astype(float)
    for tag, Xu in [("f", X),
                    ("s", X + rng.normal(0.0, MM_HIT_SMEAR_MM, X.shape))]:
        mean = np.stack([np.bincount(inv, Xu[:, k]) for k in range(3)],
                        axis=1) / n[:, None]
        Xc = Xu - mean[inv]
        M = np.empty((len(eids), 3, 3))
        for a in range(3):
            for b in range(a, 3):
                M[:, a, b] = np.bincount(inv, Xc[:, a] * Xc[:, b])
                M[:, b, a] = M[:, a, b]
        _, vecs = np.linalg.eigh(M)
        dirs = vecs[:, :, 2]                      # largest-eigenvalue axis
        sgn = np.sign((dirs * d_first).sum(axis=1))
        sgn[sgn == 0] = 1.0
        dirs = dirs * sgn[:, None]
        dirs[n < 3] = np.nan
        out[[f"{tag}_px", f"{tag}_py", f"{tag}_pz"]] = dirs
    return out


def _process_hits(df):
    """
    Convert a complete-events hits DataFrame into a per-event summary.

    Returns a DataFrame with columns:
      - em_in_<layer>, ep_in_<layer>  bool: did this particle enter the layer?
      - em_trig, ep_trig             bool: PlasticScint+LiqScint_1 on any arm
      - single_trig                  bool: >=1 arm with PlScint^LS1 (any particle)
      - double_trig                  bool: >=2 arms with PlScint^LS1 (any particle)
      - mm_double                    bool: BOTH e- and e+ hit DriftGas (geometric ideal)
      - em_edep_ls, ep_edep_ls       float: LS-only edep [MeV]
      - em_edep_all, ep_edep_all     float: edep in all scintillator layers [MeV]
                                          (PlasticScint + LS + BackScint; MM excluded)
      - em_mm_{gx,gy,gz,px,py,pz}   float: first DriftGas hit (world coords)
      - ep_mm_*                       float: same for e+
    """
    for col in ["detType", "particle"]:
        df[col] = _decode_col(df[col])

    rng = np.random.default_rng(9973 * (int(df["eventID"].iloc[0]) + 1) % (2**31))

    prim = df[df["parentID"] == 0]
    em   = prim[prim["particle"] == "e-"]
    ep   = prim[prim["particle"] == "e+"]

    all_eids = prim["eventID"].unique()
    out = pd.DataFrame({"eventID": all_eids}).set_index("eventID")

    # ── Per-particle, per-layer hit flags ─────────────────────────────────
    for prefix, p in [("em", em), ("ep", ep)]:
        for layer in SCORED_LAYERS:
            hit_evts = p.loc[p["detType"] == layer, "eventID"].unique()
            out[f"{prefix}_in_{layer}"] = False
            out.loc[out.index.isin(hit_evts), f"{prefix}_in_{layer}"] = True

        # Per-particle trigger: PlasticScint AND LiqScint_1 on the same arm
        trig = p[p["detType"].isin(TRIGGER_LAYERS)]
        if not trig.empty:
            arm_layers = trig.groupby(["eventID", "armID"])["detType"].apply(frozenset)
            trig_eids  = arm_layers[arm_layers.apply(
                lambda s: TRIGGER_LAYERS.issubset(s))
            ].index.get_level_values("eventID").unique()
            out[f"{prefix}_trig"] = out.index.isin(trig_eids)
        else:
            out[f"{prefix}_trig"] = False

        # First DriftGas hit per event — position, direction, and arm
        mm = p[p["detType"] == "DriftGas"]
        if not mm.empty:
            first_mm = mm.groupby("eventID").first()[["gx","gy","gz","px","py","pz","armID"]]
            first_mm.columns = [f"{prefix}_mm_{c}"
                                 for c in ["gx","gy","gz","px","py","pz","arm"]]
            out = out.join(first_mm)

            # Drift-track quantities (last-hit dir, centroid, PCA fits),
            # restricted to the arm of the first DriftGas hit
            arm_first = mm.groupby("eventID")["armID"].transform("first")
            trk = _drift_track_quantities(mm[mm["armID"] == arm_first], rng)
            trk.columns = [f"{prefix}_trk_{c}" for c in trk.columns]
            out = out.join(trk)

        # First trigger-scintillator hit per event (for scint-based pointing)
        ps = p[p["detType"] == "PlasticScint"]
        if not ps.empty:
            first_ps = ps.groupby("eventID").first()[["gx","gy","gz","armID"]]
            first_ps.columns = [f"{prefix}_ps_{c}" for c in ["gx","gy","gz","arm"]]
            out = out.join(first_ps)

        # Stopping detection: find the hit with minimum KE for each event.
        # The detType of that hit is the deepest scored volume the particle reached.
        # If ke_min < STOP_KE_THRESH the particle stopped there; otherwise it escaped.
        if "ke" in p.columns and not p.empty:
            idx_min = p.groupby("eventID")["ke"].idxmin().dropna().astype(int)
            if len(idx_min) > 0:
                min_rows = p.loc[idx_min.values, ["ke", "detType"]].copy()
                min_rows.index = idx_min.index  # re-index by eventID
                out[f"{prefix}_stop_ke"]  = min_rows["ke"].reindex(out.index)
                out[f"{prefix}_stop_vol"] = min_rows["detType"].reindex(out.index)

    # ── Single + double trigger (any particle, PlScint^LS1 on >=1 or >=2 arms)
    both = prim[prim["detType"].isin(TRIGGER_LAYERS)]
    if not both.empty:
        arm_layers     = both.groupby(["eventID", "armID"])["detType"].apply(frozenset)
        arm_trig       = arm_layers[arm_layers.apply(lambda s: TRIGGER_LAYERS.issubset(s))]
        arms_per_event = arm_trig.groupby("eventID").size()
        out["single_trig"] = out.index.isin(arms_per_event[arms_per_event >= 1].index)
        out["double_trig"] = out.index.isin(arms_per_event[arms_per_event >= 2].index)
    else:
        out["single_trig"] = False
        out["double_trig"] = False

    # ── Calorimetric edep: ancestry-attributed, all particles ────────────
    # Trace each hit's parentID chain back to the root primary (e- or e+).
    # This correctly attributes secondary EM shower energy — delta electrons,
    # bremsstrahlung gammas, Compton electrons, e+ annihilation gammas — to
    # the particle that produced them, regardless of detector arm.
    #
    # Geant4 assigns trackID sequentially so trackID > parentID always holds.
    # A single forward pass (sorted by trackID) propagates the root primary
    # to every descendant in O(unique tracks) with no special-casing needed.
    _NON_MM = ~df["detType"].isin(frozenset(["DriftGas", "AmpGas"]))

    # Build (eventID, trackID) → root_primary_trackID via forward pass
    _td = (df[["eventID", "trackID", "parentID"]]
           .drop_duplicates()
           .sort_values(["eventID", "trackID"]))
    _root: dict = {}
    for _r in _td.itertuples(index=False):
        _k = (_r.eventID, _r.trackID)
        if _r.parentID == 0:
            _root[_k] = _r.trackID
        else:
            _root[_k] = _root.get((_r.eventID, _r.parentID), _r.parentID)

    # Assign root primary trackID to each hit
    df["_prim_tid"] = [_root.get((e, t), t)
                       for e, t in zip(df["eventID"].values, df["trackID"].values)]

    # Map root trackID → particle type ("e-" or "e+") using primary track records
    _ptype_map = (prim[["eventID", "trackID", "particle"]]
                  .drop_duplicates()
                  .rename(columns={"trackID": "_prim_tid", "particle": "_prim_ptype"}))
    df = df.merge(_ptype_map, on=["eventID", "_prim_tid"], how="left")
    # "_prim_ptype" is "e-", "e+", or NaN for hits whose ancestry chain
    # passes through an unscored volume (rare; those hits are simply skipped)

    # Recompute _NON_MM after the merge resets the DataFrame index
    _NON_MM = ~df["detType"].isin(frozenset(["DriftGas", "AmpGas"]))

    for prefix, ptype in [("em", "e-"), ("ep", "e+")]:
        _sel = df["_prim_ptype"] == ptype
        ls_dep  = (df[_sel & df["detType"].isin(LS_LAYERS)]
                   .groupby("eventID")["edep"].sum() / 1e6)
        bs_dep  = (df[_sel & df["detType"].isin(BS_LAYERS)]
                   .groupby("eventID")["edep"].sum() / 1e6)
        all_dep = (df[_sel & _NON_MM]
                   .groupby("eventID")["edep"].sum() / 1e6)
        out[f"{prefix}_edep_ls"]  = ls_dep.reindex(out.index,  fill_value=0.0)
        out[f"{prefix}_edep_bs"]  = bs_dep.reindex(out.index,  fill_value=0.0)
        out[f"{prefix}_edep_all"] = all_dep.reindex(out.index, fill_value=0.0)

    # MM Double: both e- and e+ hit DriftGas — geometric ideal for track reconstruction
    out["mm_double"] = out["em_in_DriftGas"].eq(True) & out["ep_in_DriftGas"].eq(True)

    # MM topology: same-arm (both on one arm) vs diff-arm (on separate arms)
    _em_arm = (out["em_mm_arm"].fillna(-2)
               if "em_mm_arm" in out.columns
               else pd.Series(-2.0, index=out.index))
    _ep_arm = (out["ep_mm_arm"].fillna(-2)
               if "ep_mm_arm" in out.columns
               else pd.Series(-2.0, index=out.index))
    _both_mm = out["mm_double"]
    out["mm_same_arm"] = _both_mm & (_em_arm == _ep_arm)
    out["mm_diff_arm"]  = _both_mm & (_em_arm != _ep_arm)

    return out.reset_index()


# ── Accumulator ───────────────────────────────────────────────────────────────

class Accumulator:
    """Accumulates per-file results into histograms, never storing all rows."""

    def __init__(self):
        # Acceptance counts
        self.n_total            = {0: 0, 1: 0}
        self.n_single_trig      = {0: 0, 1: 0}
        self.n_double_trig      = {0: 0, 1: 0}
        self.n_any_mm           = {0: 0, 1: 0}
        # MM Double: both e- and e+ have DriftGas hits (geometric ideal)
        self.n_mm_double        = {0: 0, 1: 0}
        # Scintillator trigger efficiency conditioned on MM Double
        self.n_mm_double_single = {0: 0, 1: 0}
        self.n_mm_double_double = {0: 0, 1: 0}

        # MM topology: same-arm vs diff-arm subsets of mm_double
        self.n_mm_same_arm         = {0: 0, 1: 0}
        self.n_mm_diff_arm         = {0: 0, 1: 0}
        # Trigger breakdown for same-arm events (can only fire Singles)
        self.n_same_arm_single     = {0: 0, 1: 0}
        self.n_same_arm_no_trig    = {0: 0, 1: 0}
        # Trigger breakdown for diff-arm events
        self.n_diff_arm_double      = {0: 0, 1: 0}
        self.n_diff_arm_single_only = {0: 0, 1: 0}   # Singles fired but not Doubles
        self.n_diff_arm_no_trig     = {0: 0, 1: 0}

        self.layer_accept  = {0: {l: 0 for l in SCORED_LAYERS},
                              1: {l: 0 for l in SCORED_LAYERS}}

        # Invariant mass
        self.mass_bins = np.linspace(0, 25, 101)
        self.h_mass_truth        = {0: np.zeros(100), 1: np.zeros(100)}
        self.h_mass_reco         = {0: np.zeros(100), 1: np.zeros(100)}
        self.h_mass_reco_stopped = {0: np.zeros(100), 1: np.zeros(100)}
        # Reco mass using LS+BackScint energy and full-scintillator-stack energy
        self.h_mass_reco_lsbs    = {0: np.zeros(100), 1: np.zeros(100)}
        self.h_mass_reco_all     = {0: np.zeros(100), 1: np.zeros(100)}

        # Calorimeter QA: 2D histograms (true E_total = KE + mₑc², edep)
        # x-axis is E_total; 0.5 MeV bins 0-22 MeV
        _n = 44
        self.qa_ke_bins   = np.linspace(0, 22, _n + 1)   # represents E_total on x-axis
        self.qa_edep_bins = np.linspace(0, 22, _n + 1)
        self.h_qa_ls_em   = {0: np.zeros((_n, _n)), 1: np.zeros((_n, _n))}
        self.h_qa_ls_ep   = {0: np.zeros((_n, _n)), 1: np.zeros((_n, _n))}
        # LS + BackScint edep (only for particles that hit BackScint)
        self.h_qa_lsbs_em = {0: np.zeros((_n, _n)), 1: np.zeros((_n, _n))}
        self.h_qa_lsbs_ep = {0: np.zeros((_n, _n)), 1: np.zeros((_n, _n))}
        self.h_qa_all_em  = {0: np.zeros((_n, _n)), 1: np.zeros((_n, _n))}
        self.h_qa_all_ep  = {0: np.zeros((_n, _n)), 1: np.zeros((_n, _n))}
        # Stopping / BackScint counts per E_total bin
        self.h_stop_reach_em = {0: np.zeros(_n), 1: np.zeros(_n)}
        self.h_stop_stop_em  = {0: np.zeros(_n), 1: np.zeros(_n)}
        self.h_stop_reach_ep = {0: np.zeros(_n), 1: np.zeros(_n)}
        self.h_stop_stop_ep  = {0: np.zeros(_n), 1: np.zeros(_n)}
        self.h_bs_hit_em     = {0: np.zeros(_n),  1: np.zeros(_n)}
        self.h_bs_hit_ep     = {0: np.zeros(_n),  1: np.zeros(_n)}
        # Truly stopped in BS (stop_ke < STOP_KE_THRESH, stop_vol in BS_LAYERS)
        self.h_stopped_in_bs_em = {0: np.zeros(_n), 1: np.zeros(_n)}
        self.h_stopped_in_bs_ep = {0: np.zeros(_n), 1: np.zeros(_n)}
        # Escaped from LS (reached LS, not stopped in LS or BS = lateral escape)
        self.h_escaped_ls_em    = {0: np.zeros(_n), 1: np.zeros(_n)}
        self.h_escaped_ls_ep    = {0: np.zeros(_n), 1: np.zeros(_n)}

        # Upstream energy loss: E_total - edep (for particles reaching LS)
        # 2D: (E_total bin, upstream_loss bin)
        # _ul_ls  = E_total − edep_ls  (includes PlasticScint + MM + unscored material)
        # _ul_all = E_total − edep_all (only unscored material: CFRP, air, structural support)
        #           where edep_all = PlasticScint + LS1 + LS2 + BackScintL + BackScintR (no MM)
        _n_ul = _n   # same x-axis bins as calo QA (E_total)
        self.ul_loss_bins    = np.linspace(0, 18, 73)   # 0-18 MeV, 0.25 MeV/bin
        self.h_ul_ls_em      = {0: np.zeros((_n_ul, 72)), 1: np.zeros((_n_ul, 72))}
        self.h_ul_ls_ep      = {0: np.zeros((_n_ul, 72)), 1: np.zeros((_n_ul, 72))}
        # LS + PlasticScint (= edep_all − edep_bs = no BackScint, no MM)
        self.h_qa_lsplsc_em  = {0: np.zeros((_n, _n)),   1: np.zeros((_n, _n))}
        self.h_qa_lsplsc_ep  = {0: np.zeros((_n, _n)),   1: np.zeros((_n, _n))}
        self.h_ul_lsplsc_em  = {0: np.zeros((_n_ul, 72)), 1: np.zeros((_n_ul, 72))}
        self.h_ul_lsplsc_ep  = {0: np.zeros((_n_ul, 72)), 1: np.zeros((_n_ul, 72))}
        self.h_ul_all_em     = {0: np.zeros((_n_ul, 72)), 1: np.zeros((_n_ul, 72))}
        self.h_ul_all_ep     = {0: np.zeros((_n_ul, 72)), 1: np.zeros((_n_ul, 72))}

        # Energy asymmetry
        self.asym_bins        = np.linspace(0, 1, 51)
        self.h_asym_all       = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_mm        = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_mm_double = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_single    = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_double    = {0: np.zeros(50), 1: np.zeros(50)}

        # Asymmetry histograms for MM topology classes
        self.h_asym_same_arm          = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_same_arm_single   = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_same_arm_notrig   = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_diff_arm          = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_diff_arm_double   = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_diff_arm_singonly = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_diff_arm_notrig   = {0: np.zeros(50), 1: np.zeros(50)}

        # DCA of back-projected MM track to true vertex, binned by KE
        # DCA stored in cm; coarse bins for the per-KE distribution plot,
        # fine bins for the smooth median-vs-KE curve
        self.dca_bins      = np.linspace(0, 10, 101)          # 0-10 cm, 0.1 cm/bin
        self.ke_bins       = np.array([0,1,2,3,4,5,6,8,10,15])  # coarse — colours in dist plot
        self.ke_fine_bins  = np.linspace(0, 20, 81)           # 0.25 MeV/bin — median curve
        n_ke_c = len(self.ke_bins) - 1
        n_ke_f = len(self.ke_fine_bins) - 1
        self.h_dca      = {0: np.zeros((n_ke_c, 100)), 1: np.zeros((n_ke_c, 100))}
        self.h_dca_fine = {0: np.zeros((n_ke_f, 100)), 1: np.zeros((n_ke_f, 100))}

        # Opening angle resolution
        self.open_bins  = np.linspace(60, 180, 61)
        self.delta_bins = np.linspace(-60, 60, 121)
        self.h_delta_open  = {0: np.zeros((60, 120)), 1: np.zeros((60, 120))}
        self.open_rms_sum  = {0: np.zeros(60), 1: np.zeros(60)}
        self.open_rms_n    = {0: np.zeros(60), 1: np.zeros(60)}
        # 2D truth vs reco opening angle (2° bins, 0-180°)
        self.open2d_bins = np.linspace(0, 180, 91)
        self.h_open_2d   = {0: np.zeros((90, 90)), 1: np.zeros((90, 90))}

        # ── Per-particle direction-error estimators ────────────────────────
        # ψ = space angle between estimator direction and truth-at-vertex
        # (or between two estimators, see PSI_LABEL), binned vs particle KE
        n_ke_f2 = len(self.ke_fine_bins) - 1
        self.psi_bins = np.linspace(0, 45, 181)           # 0.25°/bin
        n_psi = len(self.psi_bins) - 1
        self.h_psi = {nm: {0: np.zeros((n_ke_f2, n_psi)),
                           1: np.zeros((n_ke_f2, n_psi))}
                      for nm in PSI_ESTIMATORS}

        # ── Scattering localization ────────────────────────────────────────
        # Effective scatter distance from vertex along the truth direction:
        #   s_eff = t − δ/tan(ψ)   with t = (P−V)·d̂, δ = perp displacement of
        # the first MM hit from the truth ray, ψ = upstream direction error.
        # Single-scatter exact; for distributed MS it is an effective lever arm.
        self.seff_bins = np.linspace(-150, 350, 126)      # mm, 4 mm/bin
        self.h_seff = {0: np.zeros((n_ke_f2, 125)), 1: np.zeros((n_ke_f2, 125))}

        # ── Opening-angle method comparison ────────────────────────────────
        self.mopen_bins = np.linspace(0, 180, 61)         # 3° truth bins
        n_mo, n_dl = len(self.mopen_bins) - 1, len(self.delta_bins) - 1
        self.h_mopen_delta = {mth: {0: np.zeros(n_dl), 1: np.zeros(n_dl)}
                              for mth in OPEN_METHODS}
        self.mopen_rms_sum = {mth: {0: np.zeros(n_mo), 1: np.zeros(n_mo)}
                              for mth in OPEN_METHODS}
        self.mopen_rms_n   = {mth: {0: np.zeros(n_mo), 1: np.zeros(n_mo)}
                              for mth in OPEN_METHODS}
        # Δθ vs min(KE⁻, KE⁺) — resolution gain from energy selection
        self.h_mopen_minke = {mth: {0: np.zeros((n_ke_f2, n_dl)),
                                    1: np.zeros((n_ke_f2, n_dl))}
                              for mth in OPEN_METHODS}

    def merge(self, other):
        """Add another Accumulator into this one (used after parallel processing)."""
        for et in [0, 1]:
            self.n_total[et]            += other.n_total[et]
            self.n_single_trig[et]      += other.n_single_trig[et]
            self.n_double_trig[et]      += other.n_double_trig[et]
            self.n_any_mm[et]           += other.n_any_mm[et]
            self.n_mm_double[et]        += other.n_mm_double[et]
            self.n_mm_double_single[et] += other.n_mm_double_single[et]
            self.n_mm_double_double[et] += other.n_mm_double_double[et]
            self.n_mm_same_arm[et]          += other.n_mm_same_arm[et]
            self.n_mm_diff_arm[et]          += other.n_mm_diff_arm[et]
            self.n_same_arm_single[et]      += other.n_same_arm_single[et]
            self.n_same_arm_no_trig[et]     += other.n_same_arm_no_trig[et]
            self.n_diff_arm_double[et]      += other.n_diff_arm_double[et]
            self.n_diff_arm_single_only[et] += other.n_diff_arm_single_only[et]
            self.n_diff_arm_no_trig[et]     += other.n_diff_arm_no_trig[et]
            for layer in SCORED_LAYERS:
                self.layer_accept[et][layer] += other.layer_accept[et][layer]
            for attr in [
                "h_mass_truth", "h_mass_reco", "h_mass_reco_stopped",
                "h_mass_reco_lsbs", "h_mass_reco_all",
                "h_qa_ls_em", "h_qa_ls_ep",
                "h_qa_lsbs_em", "h_qa_lsbs_ep",
                "h_qa_all_em", "h_qa_all_ep",
                "h_stop_reach_em", "h_stop_stop_em",
                "h_stop_reach_ep", "h_stop_stop_ep",
                "h_bs_hit_em", "h_bs_hit_ep",
                "h_stopped_in_bs_em", "h_stopped_in_bs_ep",
                "h_escaped_ls_em", "h_escaped_ls_ep",
                "h_ul_ls_em", "h_ul_ls_ep",
                "h_qa_lsplsc_em", "h_qa_lsplsc_ep",
                "h_ul_lsplsc_em", "h_ul_lsplsc_ep",
                "h_ul_all_em", "h_ul_all_ep",
                "h_asym_all", "h_asym_mm", "h_asym_mm_double",
                "h_asym_single", "h_asym_double",
                "h_asym_same_arm", "h_asym_same_arm_single", "h_asym_same_arm_notrig",
                "h_asym_diff_arm", "h_asym_diff_arm_double",
                "h_asym_diff_arm_singonly", "h_asym_diff_arm_notrig",
                "h_dca", "h_dca_fine", "h_delta_open", "open_rms_sum", "open_rms_n",
                "h_open_2d", "h_seff",
            ]:
                getattr(self, attr)[et] += getattr(other, attr)[et]
            for attr in ["h_psi", "h_mopen_delta", "mopen_rms_sum",
                         "mopen_rms_n", "h_mopen_minke"]:
                a, b = getattr(self, attr), getattr(other, attr)
                for key in a:
                    a[key][et] += b[key][et]

    def update(self, merged):
        """Update accumulators from a merged (EventTree + hit summary) DataFrame."""
        for et in [0, 1]:
            m = merged[merged["event_type"] == et]
            if m.empty:
                continue
            n = len(m)
            self.n_total[et] += n

            # NaN-safe bool helper: after left-merge, missing hit columns are NaN
            def _bc(col):
                return m.get(col, pd.Series(False, index=m.index)).eq(True)

            single_mask     = _bc("single_trig")
            double_mask     = _bc("double_trig")
            mm_mask         = _bc("em_in_DriftGas") | _bc("ep_in_DriftGas")
            mm_double_mask  = _bc("mm_double")   # both e- and e+ hit DriftGas
            same_arm_mask   = _bc("mm_same_arm")  # both hit MM on the SAME arm
            diff_arm_mask   = _bc("mm_diff_arm")  # each particle on a DIFFERENT arm
            single_only_m   = single_mask & ~double_mask  # single but not double

            self.n_single_trig[et]      += int(single_mask.sum())
            self.n_double_trig[et]      += int(double_mask.sum())
            self.n_any_mm[et]           += int(mm_mask.sum())
            self.n_mm_double[et]        += int(mm_double_mask.sum())
            self.n_mm_double_single[et] += int((mm_double_mask & single_mask).sum())
            self.n_mm_double_double[et] += int((mm_double_mask & double_mask).sum())

            self.n_mm_same_arm[et]          += int(same_arm_mask.sum())
            self.n_mm_diff_arm[et]          += int(diff_arm_mask.sum())
            self.n_same_arm_single[et]      += int((same_arm_mask & single_mask).sum())
            self.n_same_arm_no_trig[et]     += int((same_arm_mask & ~single_mask).sum())
            self.n_diff_arm_double[et]      += int((diff_arm_mask & double_mask).sum())
            self.n_diff_arm_single_only[et] += int((diff_arm_mask & single_only_m).sum())
            self.n_diff_arm_no_trig[et]     += int((diff_arm_mask & ~single_mask).sum())

            for layer in SCORED_LAYERS:
                self.layer_accept[et][layer] += int(
                    (_bc(f"em_in_{layer}") | _bc(f"ep_in_{layer}")).sum())

            # ── Invariant mass ────────────────────────────────────────────
            np.add.at(self.h_mass_truth[et],
                      np.clip(np.digitize(m["inv_mass"].values, self.mass_bins) - 1,
                              0, 99), 1)

            mm_cols_em = [f"em_mm_{c}" for c in ["px","py","pz"]]
            mm_cols_ep = [f"ep_mm_{c}" for c in ["px","py","pz"]]
            has_mm = all(c in m.columns for c in mm_cols_em + mm_cols_ep)

            def _reco_mass(sub):
                E_em = (sub["em_edep_ls"] + ME_MEV).values
                E_ep = (sub["ep_edep_ls"] + ME_MEV).values
                cos_th = (sub["em_mm_px"]*sub["ep_mm_px"] +
                          sub["em_mm_py"]*sub["ep_mm_py"] +
                          sub["em_mm_pz"]*sub["ep_mm_pz"]).clip(-1, 1).values
                return _inv_mass_exact(E_em, E_ep, cos_th)

            has_ls = (m.get("em_edep_ls", 0) > 0) & (m.get("ep_edep_ls", 0) > 0)
            if has_ls.any() and has_mm:
                r = m[has_ls].dropna(subset=mm_cols_em + mm_cols_ep)
                if not r.empty:
                    np.add.at(self.h_mass_reco[et],
                              np.clip(np.digitize(_reco_mass(r),
                                                  self.mass_bins) - 1, 0, 99), 1)

            # ── Stopping classification using min-KE hit ──────────────────
            # em_stop_ke / ep_stop_ke: KE at the deepest scored hit for each primary.
            # em_stop_vol / ep_stop_vol: which volume that final hit was in.
            _inf = pd.Series(np.inf, index=m.index)
            em_stop_ke   = m.get("em_stop_ke",  _inf).fillna(np.inf)
            ep_stop_ke   = m.get("ep_stop_ke",  _inf).fillna(np.inf)
            em_stop_vol  = m.get("em_stop_vol", pd.Series("", index=m.index)).fillna("")
            ep_stop_vol  = m.get("ep_stop_vol", pd.Series("", index=m.index)).fillna("")

            em_stopped_ls = (em_stop_ke < STOP_KE_THRESH) & em_stop_vol.isin(LS_LAYERS)
            ep_stopped_ls = (ep_stop_ke < STOP_KE_THRESH) & ep_stop_vol.isin(LS_LAYERS)
            em_stopped_bs = (em_stop_ke < STOP_KE_THRESH) & em_stop_vol.isin(BS_LAYERS)
            ep_stopped_bs = (ep_stop_ke < STOP_KE_THRESH) & ep_stop_vol.isin(BS_LAYERS)
            # Escaped from LS: primary reached LS but did not stop in LS or BS
            em_reached_ls = _bc("em_in_LiqScint_1") | _bc("em_in_LiqScint_2")
            ep_reached_ls = _bc("ep_in_LiqScint_1") | _bc("ep_in_LiqScint_2")
            em_escaped_ls = em_reached_ls & ~em_stopped_ls & ~em_stopped_bs
            ep_escaped_ls = ep_reached_ls & ~ep_stopped_ls & ~ep_stopped_bs

            # Keep a combined "stopped" alias used by downstream reco mass fills
            em_stopped = em_stopped_ls | em_stopped_bs
            ep_stopped = ep_stopped_ls | ep_stopped_bs
            both_stopped = em_stopped_ls & ep_stopped_ls & has_ls
            if both_stopped.any() and has_mm:
                r2 = m[both_stopped].dropna(subset=mm_cols_em + mm_cols_ep)
                if not r2.empty:
                    np.add.at(self.h_mass_reco_stopped[et],
                              np.clip(np.digitize(_reco_mass(r2),
                                                  self.mass_bins) - 1, 0, 99), 1)

            # ── LS + BackScint reco mass ───────────────────────────────────
            em_lsbs_s = (m.get("em_edep_ls", pd.Series(0.0, index=m.index)).fillna(0) +
                         m.get("em_edep_bs", pd.Series(0.0, index=m.index)).fillna(0))
            ep_lsbs_s = (m.get("ep_edep_ls", pd.Series(0.0, index=m.index)).fillna(0) +
                         m.get("ep_edep_bs", pd.Series(0.0, index=m.index)).fillna(0))
            has_lsbs  = (em_lsbs_s > 0) & (ep_lsbs_s > 0)
            if has_lsbs.any() and has_mm:
                r_l = m[has_lsbs].dropna(subset=mm_cols_em + mm_cols_ep)
                if not r_l.empty:
                    E_em = em_lsbs_s.loc[r_l.index] + ME_MEV
                    E_ep = ep_lsbs_s.loc[r_l.index] + ME_MEV
                    cos_t = (r_l["em_mm_px"]*r_l["ep_mm_px"] +
                             r_l["em_mm_py"]*r_l["ep_mm_py"] +
                             r_l["em_mm_pz"]*r_l["ep_mm_pz"]).clip(-1, 1).values
                    m_lsbs = _inv_mass_exact(E_em.values, E_ep.values, cos_t)
                    np.add.at(self.h_mass_reco_lsbs[et],
                              np.clip(np.digitize(m_lsbs, self.mass_bins) - 1, 0, 99), 1)

            # ── All scintillators reco mass ────────────────────────────────────────
            em_all_s = m.get("em_edep_all", pd.Series(0.0, index=m.index)).fillna(0)
            ep_all_s = m.get("ep_edep_all", pd.Series(0.0, index=m.index)).fillna(0)
            has_all  = (em_all_s > 0) & (ep_all_s > 0)
            if has_all.any() and has_mm:
                r_a = m[has_all].dropna(subset=mm_cols_em + mm_cols_ep)
                if not r_a.empty:
                    E_em = em_all_s.loc[r_a.index] + ME_MEV
                    E_ep = ep_all_s.loc[r_a.index] + ME_MEV
                    cos_t = (r_a["em_mm_px"]*r_a["ep_mm_px"] +
                             r_a["em_mm_py"]*r_a["ep_mm_py"] +
                             r_a["em_mm_pz"]*r_a["ep_mm_pz"]).clip(-1, 1).values
                    m_all = _inv_mass_exact(E_em.values, E_ep.values, cos_t)
                    np.add.at(self.h_mass_reco_all[et],
                              np.clip(np.digitize(m_all, self.mass_bins) - 1, 0, 99), 1)

            # ── Calorimeter QA ────────────────────────────────────────────
            for (ke_col, edep_ls_col, edep_bs_col, edep_all_col,
                 h_ls, h_lsbs, h_all, h_reach, h_stop_ls,
                 h_bs_hit, h_stopped_in_bs, h_escaped_ls,
                 h_ul_ls, h_qa_lsplsc, h_ul_lsplsc, h_ul_all,
                 is_reached_prim,
                 is_stopped_ls, is_stopped_bs, is_escaped_ls) in [
                ("em_ke", "em_edep_ls", "em_edep_bs", "em_edep_all",
                 self.h_qa_ls_em[et], self.h_qa_lsbs_em[et], self.h_qa_all_em[et],
                 self.h_stop_reach_em[et], self.h_stop_stop_em[et],
                 self.h_bs_hit_em[et],
                 self.h_stopped_in_bs_em[et], self.h_escaped_ls_em[et],
                 self.h_ul_ls_em[et],
                 self.h_qa_lsplsc_em[et], self.h_ul_lsplsc_em[et],
                 self.h_ul_all_em[et],
                 em_reached_ls,
                 em_stopped_ls, em_stopped_bs, em_escaped_ls),
                ("ep_ke", "ep_edep_ls", "ep_edep_bs", "ep_edep_all",
                 self.h_qa_ls_ep[et], self.h_qa_lsbs_ep[et], self.h_qa_all_ep[et],
                 self.h_stop_reach_ep[et], self.h_stop_stop_ep[et],
                 self.h_bs_hit_ep[et],
                 self.h_stopped_in_bs_ep[et], self.h_escaped_ls_ep[et],
                 self.h_ul_ls_ep[et],
                 self.h_qa_lsplsc_ep[et], self.h_ul_lsplsc_ep[et],
                 self.h_ul_all_ep[et],
                 ep_reached_ls,
                 ep_stopped_ls, ep_stopped_bs, ep_escaped_ls),
            ]:
                ke        = m[ke_col].values
                e_tot     = ke + ME_MEV   # total energy for x-axis
                ki        = np.clip(np.digitize(e_tot, self.qa_ke_bins) - 1,
                                    0, len(self.qa_ke_bins) - 2)
                edep_ls   = m.get(edep_ls_col,  pd.Series(0.0, index=m.index)).fillna(0).values
                edep_bs   = m.get(edep_bs_col,  pd.Series(0.0, index=m.index)).fillna(0).values
                edep_all  = m.get(edep_all_col, pd.Series(0.0, index=m.index)).fillna(0).values
                edep_lsbs = edep_ls + edep_bs

                # reached_prim: primary particle hit LS (the correct denominator for
                # stopping fractions and BackScint coverage).  A small number of events
                # have edep_ls > 0 only from secondaries while the primary stopped in
                # upstream material — using edep_ls > 0 as denominator would place those
                # in h_reach but in none of the three stopping categories, breaking sum = 1.
                reached_prim = is_reached_prim.values
                reached_any  = edep_ls > 0   # any ancestry-based LS deposit (for 2D QA)
                bs_hit       = edep_bs > 0

                if reached_any.any():
                    ei_ls  = np.clip(np.digitize(edep_ls[reached_any],  self.qa_edep_bins) - 1,
                                     0, len(self.qa_edep_bins) - 2)
                    ei_all = np.clip(np.digitize(edep_all[reached_any], self.qa_edep_bins) - 1,
                                     0, len(self.qa_edep_bins) - 2)
                    np.add.at(h_ls,  (ki[reached_any], ei_ls),  1)
                    np.add.at(h_all, (ki[reached_any], ei_all), 1)
                if reached_prim.any():
                    np.add.at(h_reach, ki[reached_prim], 1)
                if bs_hit.any():
                    ei_lsbs = np.clip(np.digitize(edep_lsbs[bs_hit], self.qa_edep_bins) - 1,
                                      0, len(self.qa_edep_bins) - 2)
                    np.add.at(h_lsbs,   (ki[bs_hit], ei_lsbs), 1)
                    np.add.at(h_bs_hit, ki[bs_hit], 1)

                # Clean stopping categories from min-KE detection (no threshold proxy)
                np.add.at(h_stop_ls,       ki[is_stopped_ls.values],  1)
                np.add.at(h_stopped_in_bs, ki[is_stopped_bs.values],  1)
                np.add.at(h_escaped_ls,    ki[is_escaped_ls.values],  1)

                # Upstream energy loss + LS+PlasticScint combined 2D QA
                # edep_lsplsc = LS1 + LS2 + PlasticScint = edep_all − edep_bs
                # Use reached_prim so the upstream loss fit only includes particles
                # that actually traversed the LS (not secondary-only deposits).
                edep_lsplsc = edep_all - edep_bs
                if reached_prim.any():
                    ul_ls     = np.clip(e_tot[reached_prim] - edep_ls[reached_prim],     0.0, None)
                    ul_lsplsc = np.clip(e_tot[reached_prim] - edep_lsplsc[reached_prim], 0.0, None)
                    ul_all    = np.clip(e_tot[reached_prim] - edep_all[reached_prim],    0.0, None)
                    n_ul      = len(self.ul_loss_bins) - 1
                    li_ls     = np.clip(np.digitize(ul_ls,     self.ul_loss_bins) - 1, 0, n_ul - 1)
                    li_lsplsc = np.clip(np.digitize(ul_lsplsc, self.ul_loss_bins) - 1, 0, n_ul - 1)
                    li_all    = np.clip(np.digitize(ul_all,    self.ul_loss_bins) - 1, 0, n_ul - 1)
                    np.add.at(h_ul_ls,    (ki[reached_prim], li_ls),     1)
                    np.add.at(h_ul_lsplsc,(ki[reached_prim], li_lsplsc), 1)
                    np.add.at(h_ul_all,   (ki[reached_prim], li_all),    1)
                    # LS+PlScint 2D QA (same axes as h_qa_ls)
                    ei_lsplsc = np.clip(
                        np.digitize(edep_lsplsc[reached_prim], self.qa_edep_bins) - 1,
                        0, len(self.qa_edep_bins) - 2)
                    np.add.at(h_qa_lsplsc, (ki[reached_prim], ei_lsplsc), 1)

            # ── Energy asymmetry ──────────────────────────────────────────
            asym = (np.abs(m["em_ke"] - m["ep_ke"]) /
                    (m["em_ke"] + m["ep_ke"]).clip(lower=1e-9))
            idx_all = np.clip(np.digitize(asym.values, self.asym_bins) - 1,
                              0, len(self.asym_bins) - 2)
            np.add.at(self.h_asym_all[et], idx_all, 1)
            if mm_mask.any():
                np.add.at(self.h_asym_mm[et],        idx_all[mm_mask.values],        1)
            if mm_double_mask.any():
                np.add.at(self.h_asym_mm_double[et], idx_all[mm_double_mask.values], 1)
            if single_mask.any():
                np.add.at(self.h_asym_single[et], idx_all[single_mask.values], 1)
            if double_mask.any():
                np.add.at(self.h_asym_double[et], idx_all[double_mask.values], 1)

            # ── MM topology asymmetry fills ───────────────────────────────
            if same_arm_mask.any():
                np.add.at(self.h_asym_same_arm[et], idx_all[same_arm_mask.values], 1)
                sa_s = (same_arm_mask & single_mask).values
                sa_n = (same_arm_mask & ~single_mask).values
                if sa_s.any():
                    np.add.at(self.h_asym_same_arm_single[et], idx_all[sa_s], 1)
                if sa_n.any():
                    np.add.at(self.h_asym_same_arm_notrig[et], idx_all[sa_n], 1)
            if diff_arm_mask.any():
                np.add.at(self.h_asym_diff_arm[et], idx_all[diff_arm_mask.values], 1)
                da_d = (diff_arm_mask & double_mask).values
                da_s = (diff_arm_mask & single_only_m).values
                da_n = (diff_arm_mask & ~single_mask).values
                if da_d.any():
                    np.add.at(self.h_asym_diff_arm_double[et],   idx_all[da_d], 1)
                if da_s.any():
                    np.add.at(self.h_asym_diff_arm_singonly[et], idx_all[da_s], 1)
                if da_n.any():
                    np.add.at(self.h_asym_diff_arm_notrig[et],   idx_all[da_n], 1)

            # ── DCA pointing ──────────────────────────────────────────────
            mm_pos_em = [f"em_mm_{c}" for c in ["gx","gy","gz"]]
            mm_dir_em = [f"em_mm_{c}" for c in ["px","py","pz"]]
            mm_pos_ep = [f"ep_mm_{c}" for c in ["gx","gy","gz"]]
            mm_dir_ep = [f"ep_mm_{c}" for c in ["px","py","pz"]]
            has_mm_both = all(c in m.columns
                              for c in mm_pos_em+mm_dir_em+mm_pos_ep+mm_dir_ep)
            if has_mm_both:
                has_any_mm = m[[mm_pos_em[0], mm_pos_ep[0]]].notna().any(axis=1)
                q4 = m[has_any_mm].copy()
                if not q4.empty:
                    for pos_cols, dir_cols, ke_col in [
                        (mm_pos_em, mm_dir_em, "em_ke"),
                        (mm_pos_ep, mm_dir_ep, "ep_ke"),
                    ]:
                        if not all(c in q4.columns for c in pos_cols + dir_cols):
                            continue
                        valid = q4[pos_cols[0]].notna()
                        qv = q4[valid]
                        if qv.empty:
                            continue
                        P = qv[pos_cols].values
                        d = qv[dir_cols].values
                        V = qv[["vtx_x","vtx_y","vtx_z"]].values
                        VP  = P - V
                        t   = (VP * d).sum(axis=1, keepdims=True)
                        dca = np.linalg.norm(VP - t * d, axis=1)
                        ke_this  = qv[ke_col].values
                        dca_cm   = dca / 10.0  # mm → cm
                        dca_idx  = np.clip(np.digitize(dca_cm, self.dca_bins) - 1,
                                           0, len(self.dca_bins) - 2)
                        ke_idx_c = np.clip(np.digitize(ke_this, self.ke_bins) - 1,
                                           0, len(self.ke_bins) - 2)
                        ke_idx_f = np.clip(np.digitize(ke_this, self.ke_fine_bins) - 1,
                                           0, len(self.ke_fine_bins) - 2)
                        np.add.at(self.h_dca[et],      (ke_idx_c, dca_idx), 1)
                        np.add.at(self.h_dca_fine[et], (ke_idx_f, dca_idx), 1)

            # ── Opening angle ─────────────────────────────────────────────
            has_mm_q5 = has_mm_both and all(c in m.columns
                                             for c in mm_dir_em + mm_dir_ep)
            if has_mm_q5:
                q5 = m[m[mm_dir_em[0]].notna() & m[mm_dir_ep[0]].notna()].copy()
                if not q5.empty:
                    d_em = q5[mm_dir_em].values
                    d_ep = q5[mm_dir_ep].values
                    cos_reco    = (d_em * d_ep).sum(axis=1).clip(-1, 1)
                    theta_reco  = np.degrees(np.arccos(cos_reco))
                    theta_truth = q5["openingAngle"].values
                    delta       = theta_reco - theta_truth
                    truth_idx   = np.clip(np.digitize(theta_truth, self.open_bins) - 1,
                                          0, len(self.open_bins) - 2)
                    delta_idx   = np.clip(np.digitize(delta, self.delta_bins) - 1,
                                          0, len(self.delta_bins) - 2)
                    np.add.at(self.h_delta_open[et], (truth_idx, delta_idx), 1)
                    np.add.at(self.open_rms_sum[et], truth_idx, delta**2)
                    np.add.at(self.open_rms_n[et],   truth_idx, 1)
                    # 2D truth vs reco
                    ti2 = np.clip(np.digitize(theta_truth, self.open2d_bins) - 1, 0, 89)
                    ri2 = np.clip(np.digitize(theta_reco,  self.open2d_bins) - 1, 0, 89)
                    np.add.at(self.h_open_2d[et], (ti2, ri2), 1)

            # ── Direction-error decomposition + method comparison ─────────
            rng_u = np.random.default_rng(
                104729 * (int(m["eventID"].iloc[0]) + 1 + et) % (2**31))
            V = m[["vtx_x", "vtx_y", "vtx_z"]].values.astype(float)
            est_dirs = {}   # est_dirs[prefix][method] → (N,3) direction array
            for prefix, ke_col in [("em", "em_ke"), ("ep", "ep_ke")]:
                d_tr  = _arr(m, [f"{prefix}_{c}"       for c in ("px", "py", "pz")])
                d_f   = _arr(m, [f"{prefix}_mm_{c}"    for c in ("px", "py", "pz")])
                P     = _arr(m, [f"{prefix}_mm_{c}"    for c in ("gx", "gy", "gz")])
                d_l   = _arr(m, [f"{prefix}_trk_l_{c}" for c in ("px", "py", "pz")])
                d_pca = _arr(m, [f"{prefix}_trk_f_{c}" for c in ("px", "py", "pz")])
                d_pcs = _arr(m, [f"{prefix}_trk_s_{c}" for c in ("px", "py", "pz")])
                C     = _arr(m, [f"{prefix}_trk_c_{c}" for c in ("gx", "gy", "gz")])
                PS    = _arr(m, [f"{prefix}_ps_{c}"    for c in ("gx", "gy", "gz")])
                ps_arm = (m.get(f"{prefix}_ps_arm", pd.Series(np.nan, index=m.index))
                          .fillna(-1).values.astype(int))
                ke = m[ke_col].values
                ki = np.clip(np.digitize(ke, self.ke_fine_bins) - 1,
                             0, len(self.ke_fine_bins) - 2)

                # Smear trigger-scint hit in the two transverse directions of
                # its arm (arms 0,1: normal=X → smear y,z; arms 2,3: smear x,y)
                PS_sm = PS.copy()
                noise = rng_u.normal(0.0, SCINT_SMEAR_MM, PS.shape)
                on_x = (ps_arm == 0) | (ps_arm == 1)
                on_z = (ps_arm == 2) | (ps_arm == 3)
                PS_sm[on_x, 1] += noise[on_x, 1]; PS_sm[on_x, 2] += noise[on_x, 2]
                PS_sm[on_z, 0] += noise[on_z, 0]; PS_sm[on_z, 1] += noise[on_z, 1]

                est_dirs[prefix] = {
                    "first":   d_f,
                    "fit":     d_pca,
                    "fits":    d_pcs,
                    "vline":   P - V,
                    "vcen":    C - V,
                    "nomline": P,
                    "scint":   PS_sm - V,
                }

                # ψ fills (per-particle direction errors)
                for name, A, B in [
                    ("first",        d_tr, d_f),
                    ("last",         d_tr, d_l),
                    ("drift",        d_f,  d_l),
                    ("fit",          d_tr, d_pca),
                    ("fit_vs_first", d_f,  d_pca),
                    ("vline",        d_tr, P - V),
                    ("vcen",         d_tr, C - V),
                    ("nomline",      d_tr, P),
                    ("scint",        d_tr, PS_sm - V),
                ]:
                    psi = _angle_deg(A, B)
                    ok  = np.isfinite(psi)
                    if not ok.any():
                        continue
                    pi = np.clip(np.digitize(psi[ok], self.psi_bins) - 1,
                                 0, len(self.psi_bins) - 2)
                    np.add.at(self.h_psi[name][et], (ki[ok], pi), 1)

                # Scattering localization from first-MM-hit kinematics
                psi_up = _angle_deg(d_tr, d_f)
                t_par  = ((P - V) * d_tr).sum(axis=1)
                d_perp = np.linalg.norm((P - V) - t_par[:, None] * d_tr, axis=1)
                ok = np.isfinite(psi_up) & (psi_up > 1.5)
                if ok.any():
                    s_eff = t_par[ok] - d_perp[ok] / np.tan(np.radians(psi_up[ok]))
                    # drop under/overflows rather than piling them in edge bins
                    inr = (s_eff >= self.seff_bins[0]) & (s_eff < self.seff_bins[-1])
                    si = np.digitize(s_eff[inr], self.seff_bins) - 1
                    np.add.at(self.h_seff[et], (ki[ok][inr], si), 1)

            # Opening-angle reconstruction methods
            theta_truth_all = m["openingAngle"].values
            min_ke = np.minimum(m["em_ke"].values, m["ep_ke"].values)
            kmi = np.clip(np.digitize(min_ke, self.ke_fine_bins) - 1,
                          0, len(self.ke_fine_bins) - 2)
            for mth in OPEN_METHODS:
                A, B = est_dirs["em"][mth], est_dirs["ep"][mth]
                th_reco = _angle_deg(A, B)
                ok = np.isfinite(th_reco)
                if not ok.any():
                    continue
                delta = th_reco[ok] - theta_truth_all[ok]
                di = np.clip(np.digitize(delta, self.delta_bins) - 1,
                             0, len(self.delta_bins) - 2)
                np.add.at(self.h_mopen_delta[mth][et], di, 1)
                ti = np.clip(np.digitize(theta_truth_all[ok], self.mopen_bins) - 1,
                             0, len(self.mopen_bins) - 2)
                np.add.at(self.mopen_rms_sum[mth][et], ti, delta**2)
                np.add.at(self.mopen_rms_n[mth][et],   ti, 1)
                np.add.at(self.h_mopen_minke[mth][et], (kmi[ok], di), 1)


# ── Per-file processing ───────────────────────────────────────────────────────

def process_file(filepath, chunk_size=300_000):
    """Process one ROOT file; returns a filled Accumulator (safe for multiprocessing)."""
    accum = Accumulator()
    with uproot.open(filepath) as f:
        if "EventTree" not in f or "HitTree" not in f:
            print(f"  WARNING: missing trees in {filepath}", file=sys.stderr)
            return 0

        evt_df = f["EventTree"].arrays(library="pd")
        if evt_df.empty:
            return 0

        leftover   = None
        hit_chunks = []

        for chunk in f["HitTree"].iterate(HIT_COLS, step_size=chunk_size, library="pd"):
            # uproot returns Char[32]/C branches as AwkwardExtensionArray.
            # pd.concat([awk, awk]) falls back to element-wise Python iteration
            # and takes ~90 s per chunk.  Convert to plain numpy object arrays
            # immediately so all subsequent concats and comparisons are O(ms).
            chunk["detType"]  = np.asarray(chunk["detType"])
            chunk["particle"] = np.asarray(chunk["particle"])

            if leftover is not None:
                chunk = pd.concat([leftover, chunk], ignore_index=True)
            last_eid = chunk["eventID"].iloc[-1]
            complete = chunk[chunk["eventID"] < last_eid]
            leftover = chunk[chunk["eventID"] >= last_eid]
            if not complete.empty:
                hit_chunks.append(_process_hits(complete.copy()))
            del chunk

        if leftover is not None and not leftover.empty:
            hit_chunks.append(_process_hits(leftover))

        if not hit_chunks:
            return len(evt_df)

        hit_summary = pd.concat(hit_chunks, ignore_index=True)
        del hit_chunks

        merged = evt_df.merge(hit_summary, on="eventID", how="left")
        del evt_df, hit_summary

        accum.update(merged)
        del merged
        gc.collect()
    return accum


# ── Plotting ──────────────────────────────────────────────────────────────────

def _legend_patches():
    return [mpatches.Patch(color=ETYPE_COLOR[e], label=ETYPE_LABEL[e])
            for e in [0, 1]]


def plot_acceptance(pdf, acc):
    fig, axes = plt.subplots(1, 3, figsize=(19, 5))

    # ── Panel 0: layer acceptance ─────────────────────────────────────────
    ax = axes[0]
    layers = SCORED_LAYERS
    x = np.arange(len(layers))
    w = 0.35
    for i, et in enumerate([0, 1]):
        n = acc.n_total[et]
        if n == 0:
            continue
        fracs = [acc.layer_accept[et][l] / n for l in layers]
        ax.bar(x + (i - 0.5) * w, fracs, w, color=ETYPE_COLOR[et],
               alpha=0.85, label=ETYPE_LABEL[et], edgecolor="k", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([LAYER_LABELS[l] for l in layers], rotation=25, ha="right")
    ax.set_ylabel("Fraction of events with ≥1 hit (e⁻ or e⁺)")
    ax.set_title("Layer acceptance")
    ax.set_ylim(0, 1.1); ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.legend(); ax.grid(True, axis="y", alpha=0.3)

    # ── Panel 1: trigger rates as fraction of all events ──────────────────
    ax = axes[1]
    trig_labels = [
        "MM Double\n(both particles\nin DriftGas)",
        "Single trigger\n(PlScint∧LS1,\n≥1 arm)",
        "Double trigger\n(PlScint∧LS1,\n≥2 arms)",
    ]
    x3 = np.arange(3)
    for i, et in enumerate([0, 1]):
        n = acc.n_total[et]
        if n == 0:
            continue
        vals = [acc.n_mm_double[et]   / n,
                acc.n_single_trig[et] / n,
                acc.n_double_trig[et] / n]
        ax.bar(x3 + (i - 0.5) * 0.35, vals, 0.35, color=ETYPE_COLOR[et],
               alpha=0.85, label=ETYPE_LABEL[et], edgecolor="k", linewidth=0.4)
        for xi, vi in zip(x3 + (i - 0.5) * 0.35, vals):
            ax.text(xi, vi + 0.01, f"{vi*100:.1f}%", ha="center",
                    va="bottom", fontsize=8)
    ax.set_xticks(x3); ax.set_xticklabels(trig_labels, fontsize=9)
    ax.set_ylabel("Fraction of all events")
    ax.set_title("Trigger rates\n(fraction of generated events)")
    ax.set_ylim(0, 1.1); ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.legend(); ax.grid(True, axis="y", alpha=0.3)

    # ── Panel 2: scintillator efficiency conditioned on MM Double ─────────
    ax = axes[2]
    cond_labels = [
        "Single trigger\n(PlScint∧LS1,\n≥1 arm)",
        "Double trigger\n(PlScint∧LS1,\n≥2 arms)",
    ]
    x2 = np.arange(2)
    for i, et in enumerate([0, 1]):
        n_mm = acc.n_mm_double[et]
        if n_mm == 0:
            continue
        vals = [acc.n_mm_double_single[et] / n_mm,
                acc.n_mm_double_double[et] / n_mm]
        ax.bar(x2 + (i - 0.5) * 0.35, vals, 0.35, color=ETYPE_COLOR[et],
               alpha=0.85, label=f"{ETYPE_LABEL[et]}  (N_MM={n_mm:,})",
               edgecolor="k", linewidth=0.4)
        for xi, vi in zip(x2 + (i - 0.5) * 0.35, vals):
            ax.text(xi, vi + 0.01, f"{vi*100:.1f}%", ha="center",
                    va="bottom", fontsize=8)
    ax.set_xticks(x2); ax.set_xticklabels(cond_labels, fontsize=9)
    ax.set_ylabel("Fraction of MM-double events")
    ax.set_title("Scintillator trigger efficiency\ngiven MM Double (both tracks present)")
    ax.set_ylim(0, 1.1); ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.legend(fontsize=8); ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle("Detection Acceptance", fontsize=13)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_invariant_mass(pdf, acc):
    bins = acc.mass_bins
    cen  = 0.5 * (bins[:-1] + bins[1:])

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    ax = axes[0]
    for et in [0, 1]:
        h = acc.h_mass_truth[et]
        norm = h.sum()
        if norm == 0: continue
        ax.step(cen, h / max(norm, 1), where="mid",
                color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                label=f"{ETYPE_LABEL[et]}  (N={int(norm):,})")
    ax.axvline(16.8, color="red", lw=1.2, ls=":", label="m_X17 = 16.8 MeV")
    ax.set_xlabel("True e⁺e⁻ invariant mass  [MeV]")
    ax.set_ylabel("Normalised counts")
    ax.set_title("Truth invariant mass")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.set_ylim(0, _mass_auto_ylim(acc.h_mass_truth[0], acc.h_mass_truth[1]))

    ax = axes[1]
    for et in [0, 1]:
        h = acc.h_mass_reco[et]
        norm = h.sum()
        if norm == 0: continue
        ax.step(cen, h / norm, where="mid",
                color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                label=f"{ETYPE_LABEL[et]}  (N={int(norm):,})")
    ax.axvline(16.8, color="red", lw=1.2, ls=":", label="m_X17 = 16.8 MeV")
    ax.set_xlabel("Reco invariant mass  [MeV]")
    ax.set_title("Reco mass — LS only (uncorrected energy)\n"
                 r"exact: $M = \sqrt{2m_e^2 + 2(E_{e^-}E_{e^+} - p_{e^-}p_{e^+}\cos\theta)}$")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.set_ylim(0, _mass_auto_ylim(acc.h_mass_reco[0], acc.h_mass_reco[1]))

    ax = axes[2]
    for et in [0, 1]:
        h = acc.h_mass_reco_all[et]
        norm = h.sum()
        if norm == 0: continue
        ax.step(cen, h / norm, where="mid",
                color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                label=f"{ETYPE_LABEL[et]}  (N={int(norm):,})")
    ax.axvline(16.8, color="red", lw=1.2, ls=":", label="m_X17 = 16.8 MeV")
    ax.set_xlabel("Reco invariant mass  [MeV]")
    ax.set_title("Reco mass — all scintillators\n(PlScint + LS + BackScint energy)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.set_ylim(0, _mass_auto_ylim(acc.h_mass_reco_all[0], acc.h_mass_reco_all[1]))

    fig.suptitle("Invariant Mass  (exact formula including mₑ)", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)

    # ── Page 2: energy-corrected mass (X17 signal only) ───────────────────
    # The correction E_est = (edep_ls + b)/(1-a) requires fit coefficients from
    # the upstream-loss histograms.  These are only available AFTER all files are
    # accumulated, so we compute the corrected mass here in the plot function.
    # We apply the exact mass formula with the corrected energies.
    l_cen2 = 0.5 * (acc.ul_loss_bins[:-1] + acc.ul_loss_bins[1:])
    e_cen2 = 0.5 * (acc.qa_ke_bins[:-1]   + acc.qa_ke_bins[1:])

    def _get_corr(h_ul):
        med   = _compute_medians(h_ul, l_cen2)
        valid = (~np.isnan(med)) & (e_cen2 > 0.5)
        if valid.sum() < 4: return None
        try: return tuple(np.polyfit(e_cen2[valid], med[valid], 1))
        except: return None

    fig2, ax2 = plt.subplots(figsize=(10, 5))
    et = 0  # X17 signal
    for h, label, color, ls_style in [
        (acc.h_mass_truth[et],    "Truth",                    "#2ca02c", "-"),
        (acc.h_mass_reco[et],     "Reco — LS only",           "#e84040", "--"),
        (acc.h_mass_reco_all[et], "Reco — all scintillators", "#9467bd", "--"),
    ]:
        n = h.sum()
        if n < 1: continue
        ax2.step(cen, h / n, where="mid", color=color, lw=1.5, ls=ls_style, label=label)

    # Corrected mass: apply E_est = (edep + b)/(1-a) to mass axis.
    # Derivation: for symmetric decays E_em ≈ E_ep ≈ E_meas,
    #   M_corr ≈ M_reco * sqrt(E_em_corr * E_ep_corr) / sqrt(E_em_meas * E_ep_meas)
    #          ≈ M_reco * 1/sqrt((1-a_em)*(1-a_ep))  (large-E limit, b negligible)
    # More precisely, scale = median(E_corr/E_meas) evaluated at the peak energy.
    c_em = _get_corr(acc.h_ul_ls_em[et])
    c_ep = _get_corr(acc.h_ul_ls_ep[et])
    c_em_all = _get_corr(acc.h_ul_all_em[et])
    c_ep_all = _get_corr(acc.h_ul_all_ep[et])

    for coeffs_em, coeffs_ep, h_base, color, label in [
        (c_em,     c_ep,     acc.h_mass_reco[et],      "#1f77b4", "LS corrected"),
        (c_em_all, c_ep_all, acc.h_mass_reco_all[et],  "#7f7f7f", "All scint. corrected"),
    ]:
        if coeffs_em is None or coeffs_ep is None: continue
        a_em, _ = coeffs_em;  a_ep, _ = coeffs_ep
        if abs(1-a_em) < 0.01 or abs(1-a_ep) < 0.01: continue
        scale = np.sqrt(1.0 / ((1-a_em)*(1-a_ep)))
        n = h_base.sum()
        if n < 1: continue
        ax2.step(cen * scale, h_base / n, where="mid",
                 color=color, lw=2.5, ls="-", label=f"{label}  (×{scale:.2f})")

    ax2.axvline(16.8, color="red", lw=1.4, ls=":", label="m_X17 = 16.8 MeV")
    ax2.set_xlabel("Invariant mass  [MeV]")
    ax2.set_ylabel("Normalised counts")
    ax2.set_ylim(0, _mass_auto_ylim(acc.h_mass_reco[et], acc.h_mass_reco_all[et]))
    ax2.set_title(
        "X17 signal: energy-corrected invariant mass\n"
        r"Correction: $E_{est} = (E_{dep} + b)/(1-a)$ from upstream-loss fit  "
        r"⟹  $M_{corr} \approx M_{reco}/\sqrt{(1-a_{e^-})(1-a_{e^+})}$  [large-E approx]")
    ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3)
    fig2.tight_layout(); pdf.savefig(fig2); plt.close(fig2)


def _calo_qa_page(pdf, acc, et_sel, title_suffix):
    """One calorimeter QA page. et_sel=None → combined, 0=X17, 1=IPC."""
    ke_cen   = 0.5 * (acc.qa_ke_bins[:-1]   + acc.qa_ke_bins[1:])
    edep_cen = 0.5 * (acc.qa_edep_bins[:-1] + acc.qa_edep_bins[1:])

    def _get(d):
        """Sum over requested event types."""
        if et_sel is None:
            return sum(d.values())
        return d[et_sel].copy()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Calorimeter QA: LS energy containment  — {title_suffix}", fontsize=12)

    for col, (particle, h_ls) in enumerate([
        ("e⁻", _get(acc.h_qa_ls_em)),
        ("e⁺", _get(acc.h_qa_ls_ep)),
    ]):
        ax = axes[0, col]
        row_sums = h_ls.sum(axis=1, keepdims=True).clip(1)
        h_norm = h_ls / row_sums
        im = ax.pcolormesh(ke_cen, edep_cen, h_norm.T, cmap="viridis", vmin=0)
        ax.plot([0, 22], [0, 22], "w--", lw=1, label="perfect containment")
        medians = []
        for ki in range(len(ke_cen)):
            row = h_ls[ki]; n = row.sum()
            if n < 5: medians.append(np.nan); continue
            cdf = np.cumsum(row) / n
            medians.append(edep_cen[np.searchsorted(cdf, 0.5)])
        valid = ~np.isnan(medians)
        ax.plot(ke_cen[valid], np.array(medians)[valid], "r-o", ms=3, lw=1.5,
                label="median edep")
        plt.colorbar(im, ax=ax, label="P(LS edep | true KE)")
        ax.set_xlabel("True E_total  [MeV]"); ax.set_ylabel("LS edep  [MeV]")
        ax.set_title(f"{particle}: LS edep vs true E_total  (LiqScint_1 + LiqScint_2 only)")
        ax.legend(fontsize=8)

    ax = axes[1, 0]
    h_tot = _get(acc.h_qa_all_em) + _get(acc.h_qa_all_ep)
    row_sums = h_tot.sum(axis=1, keepdims=True).clip(1)
    h_norm = h_tot / row_sums
    im = ax.pcolormesh(ke_cen, edep_cen, h_norm.T, cmap="viridis", vmin=0)
    ax.plot([0, 22], [0, 22], "w--", lw=1, label="perfect containment")
    medians = []
    for ki in range(len(ke_cen)):
        row = h_tot[ki]; n = row.sum()
        if n < 5: medians.append(np.nan); continue
        cdf = np.cumsum(row) / n
        medians.append(edep_cen[np.searchsorted(cdf, 0.5)])
    valid = ~np.isnan(medians)
    ax.plot(ke_cen[valid], np.array(medians)[valid], "r-o", ms=3, lw=1.5, label="median")
    plt.colorbar(im, ax=ax, label="P(total edep | true KE)")
    ax.set_xlabel("True E_total  [MeV]"); ax.set_ylabel("Total scored edep  [MeV]")
    ax.set_title("e⁻ + e⁺: total edep in all scintillator layers vs true E_total")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    for particle, h_reach, h_stop, color in [
        ("e⁻", _get(acc.h_stop_reach_em), _get(acc.h_stop_stop_em), "#4a90d9"),
        ("e⁺", _get(acc.h_stop_reach_ep), _get(acc.h_stop_stop_ep), "#e84040"),
    ]:
        safe = h_reach > 5
        frac = np.where(safe, h_stop / h_reach.clip(1), np.nan)
        ax.plot(ke_cen[safe], frac[safe], color=color, lw=2, marker="o", ms=3,
                label=f"{particle} (stopped / reached LS)")
    ax.set_xlabel("True E_total  [MeV]")
    ax.set_ylabel("Fraction stopped in LS\n(reached LS, no BackScint hit)")
    ax.set_title("Containment fraction vs KE")
    ax.set_ylim(0, 1.05); ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.legend(); ax.grid(True, alpha=0.3)

    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_calorimeter_qa(pdf, acc):
    """Three QA pages: combined, X17-only, IPC-only."""
    _calo_qa_page(pdf, acc, None, "X17 + IPC combined")
    _calo_qa_page(pdf, acc, 0,    "X17 signal only")
    _calo_qa_page(pdf, acc, 1,    "IPC background only")


def _compute_medians(h2d, cen):
    """Return median of each row of a 2D histogram; NaN for rows with <5 counts."""
    medians = np.full(h2d.shape[0], np.nan)
    for i, row in enumerate(h2d):
        n = row.sum()
        if n < 5:
            continue
        cdf = np.cumsum(row) / n
        medians[i] = cen[np.searchsorted(cdf, 0.5)]
    return medians


def _mass_auto_ylim(*hists, margin=1.35):
    """Y upper limit for normalised mass plots that clips delta-function truth spikes.
    Any bin with normalised height > 0.5 is a truth spike and is excluded from the
    maximum, so the scale is set by the reco distributions.
    """
    top = 0.0
    for h in hists:
        n = h.sum()
        if n < 1:
            continue
        normed = h / n
        safe = normed[normed < 0.5]
        if len(safe):
            top = max(top, float(safe.max()))
    return top * margin if top > 0 else 0.10


def _compute_percentile(h2d, cen, p):
    """Return the p-th percentile (0-1) for each row; NaN for rows with <5 counts."""
    result = np.full(h2d.shape[0], np.nan)
    for i, row in enumerate(h2d):
        n = row.sum()
        if n < 5:
            continue
        cdf = np.cumsum(row) / n
        idx = np.searchsorted(cdf, p)
        result[i] = cen[min(idx, len(cen) - 1)]
    return result


def plot_upstream_correction(pdf, acc):
    """Upstream energy loss analysis and linear correction for LS calorimetry.

    Page 1 (2×3): Upstream-loss distributions for three detector combinations
                  (LS, LS+PlScint, all scintillators) for e⁻ and e⁺.
    Page 2 (2×3): Corrected E_estimated vs E_total for the three methods.
    Page 3 (1 wide): Invariant mass — uncorrected + all three corrected curves.

    Correction formula
    ------------------
    Fit:       upstream_loss = a · E_total + b
    Inverted:  E_total_est  = (edep + b) / (1 − a)

    The mass comparison uses the large-E approximation
      M_corr ≈ M_reco · sqrt[1 / ((1−a_em)·(1−a_ep))]
    i.e., both particles' energies are scaled by 1/(1−a).  This is exact at
    the mass peak for symmetric (equal-energy) decays, and approximate in the
    tails where the two particles have very different energies.
    """
    e_cen    = 0.5 * (acc.qa_ke_bins[:-1]   + acc.qa_ke_bins[1:])
    edep_cen = 0.5 * (acc.qa_edep_bins[:-1] + acc.qa_edep_bins[1:])
    l_cen    = 0.5 * (acc.ul_loss_bins[:-1] + acc.ul_loss_bins[1:])

    def _fit_loss(h_ul):
        """Linear fit: upstream_loss = a·E_total + b  →  returns (a, b) or None."""
        med   = _compute_medians(h_ul, l_cen)
        valid = (~np.isnan(med)) & (e_cen > 0.5)
        if valid.sum() < 4:
            return None
        try:
            return tuple(np.polyfit(e_cen[valid], med[valid], 1))
        except Exception:
            return None

    def _draw_ul_panel(ax, h_ul, particle, lbl, key, fits_dict):
        rs = h_ul.sum(axis=1, keepdims=True).clip(1)
        im = ax.pcolormesh(e_cen, l_cen, (h_ul / rs).T, cmap="plasma", vmin=0)
        plt.colorbar(im, ax=ax, label="P(loss | E_total)")
        med = _compute_medians(h_ul, l_cen)
        p16 = _compute_percentile(h_ul, l_cen, 0.16)
        p84 = _compute_percentile(h_ul, l_cen, 0.84)
        v   = ~np.isnan(med)
        ax.plot(e_cen[v], med[v], "w-", lw=2, label="Median")
        ax.fill_between(e_cen[v], p16[v], p84[v], alpha=0.25, color="white",
                        label="16–84% band")
        coeffs = _fit_loss(h_ul)
        if coeffs is not None:
            fits_dict[key] = coeffs
            a, b = coeffs
            ax.plot(e_cen, np.polyval(coeffs, e_cen), "r--", lw=2,
                    label=f"Fit: {a:.3f}·E {b:+.2f} MeV")
        ax.set_xlabel("True E_total  [MeV]")
        ax.set_ylabel("Upstream loss  [MeV]")
        ax.set_title(f"{particle} — {lbl}\nloss = E_total − edep")
        ax.set_ylim(0, acc.ul_loss_bins[-1])
        ax.legend(fontsize=8)

    def _draw_corr_panel(ax, h_qa, key, fits_dict, particle, lbl):
        rs  = h_qa.sum(axis=1, keepdims=True).clip(1)
        med = _compute_medians(h_qa, edep_cen)
        v   = ~np.isnan(med)
        if key not in fits_dict:
            ax.set_title(f"{particle} — {lbl}\n(no fit available)")
            return
        a, b = fits_dict[key]
        if abs(1 - a) < 0.01:
            ax.set_title(f"{particle} — {lbl}\nfit slope ≈ 1, correction degenerate")
            return
        corr_bins = (acc.qa_edep_bins + b) / (1 - a)
        im = ax.pcolormesh(e_cen, corr_bins[:-1], (h_qa / rs).T, cmap="viridis", vmin=0)
        ax.plot([0, 22], [0, 22], "w--", lw=1, label="perfect (y = x)")
        ax.plot(e_cen[v], (med[v] + b) / (1 - a), "r-o", ms=3, lw=1.5,
                label="corrected median")
        plt.colorbar(im, ax=ax)
        ax.set_xlabel("True E_total  [MeV]")
        ax.set_ylabel("E_est = (edep + b)/(1−a)  [MeV]")
        ax.set_title(f"{particle} — {lbl} corrected\n"
                     f"E_est = (edep {b:+.2f}) / (1−{a:.3f})")
        ax.set_ylim(0, 25)
        ax.legend(fontsize=8)

    # ── Page 1: upstream loss distributions ──────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(19, 11))
    fig.suptitle(
        "Upstream Energy Loss (X17 signal)  —  loss = E_total − deposited energy\n"
        "columns: LS only  |  LS + PlasticScint  |  all scintillators (PlScint+LS+BS)",
        fontsize=12)
    fits = {}
    for r, c, h_ul, particle, lbl, key in [
        (0, 0, acc.h_ul_ls_em[0],     "e⁻", "LS only",           "ls_em"),
        (0, 1, acc.h_ul_lsplsc_em[0], "e⁻", "LS + PlasticScint", "lsplsc_em"),
        (0, 2, acc.h_ul_all_em[0],    "e⁻", "All scintillators", "all_em"),
        (1, 0, acc.h_ul_ls_ep[0],     "e⁺", "LS only",           "ls_ep"),
        (1, 1, acc.h_ul_lsplsc_ep[0], "e⁺", "LS + PlasticScint", "lsplsc_ep"),
        (1, 2, acc.h_ul_all_ep[0],    "e⁺", "All scintillators", "all_ep"),
    ]:
        _draw_ul_panel(axes[r, c], h_ul, particle, lbl, key, fits)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    # ── Page 2: corrected E_estimated vs E_total ─────────────────────────
    fig2, axes2 = plt.subplots(2, 3, figsize=(19, 11))
    fig2.suptitle(
        r"Corrected E estimate vs True E_total  (X17 signal)   E_est = (edep + b) / (1−a)"
        "\ncolumns: LS  |  LS + PlasticScint  |  all scintillators",
        fontsize=12)
    for r, c, h_qa, key, particle, lbl in [
        (0, 0, acc.h_qa_ls_em[0],     "ls_em",     "e⁻", "LS only"),
        (0, 1, acc.h_qa_lsplsc_em[0], "lsplsc_em", "e⁻", "LS + PlasticScint"),
        (0, 2, acc.h_qa_all_em[0],    "all_em",    "e⁻", "All scintillators"),
        (1, 0, acc.h_qa_ls_ep[0],     "ls_ep",     "e⁺", "LS only"),
        (1, 1, acc.h_qa_lsplsc_ep[0], "lsplsc_ep", "e⁺", "LS + PlasticScint"),
        (1, 2, acc.h_qa_all_ep[0],    "all_ep",    "e⁺", "All scintillators"),
    ]:
        _draw_corr_panel(axes2[r, c], h_qa, key, fits, particle, lbl)
    fig2.tight_layout(); pdf.savefig(fig2); plt.close(fig2)

    # ── Page 3: corrected invariant mass comparison ───────────────────────
    mass_bins = acc.mass_bins
    mass_cen  = 0.5 * (mass_bins[:-1] + mass_bins[1:])
    et = 0  # X17 signal

    fig3, ax3 = plt.subplots(figsize=(13, 5.5))

    # Uncorrected (thin dashed)
    for h, label, color in [
        (acc.h_mass_truth[et],     "Truth",                    "#2ca02c"),
        (acc.h_mass_reco[et],      "Reco — LS only",           "#e84040"),
        (acc.h_mass_reco_lsbs[et], "Reco — LS + BackScint",    "#ff7f0e"),
        (acc.h_mass_reco_all[et],  "Reco — all scintillators", "#9467bd"),
    ]:
        n = h.sum()
        if n < 1: continue
        ax3.step(mass_cen, h / n, where="mid", color=color, lw=1.5, ls="--",
                 alpha=0.7, label=label)

    # Corrected (thick solid): rescale mass axis by sqrt[1/((1-a_em)*(1-a_ep))]
    # This is the large-E approximation: E_corr/E_meas ≈ 1/(1-a).
    for key_em, key_ep, h_base, color, label in [
        ("ls_em",     "ls_ep",     acc.h_mass_reco[et],      "#1f77b4",
         "LS corrected"),
        ("lsplsc_em", "lsplsc_ep", acc.h_mass_reco[et],      "#17becf",
         "LS+PlScint corrected"),
        ("all_em",    "all_ep",    acc.h_mass_reco_all[et],  "#7f7f7f",
         "All scintillators corrected"),
    ]:
        if key_em not in fits or key_ep not in fits: continue
        a_em, _ = fits[key_em];  a_ep, _ = fits[key_ep]
        if abs(1 - a_em) < 0.01 or abs(1 - a_ep) < 0.01: continue
        scale = np.sqrt(1.0 / ((1 - a_em) * (1 - a_ep)))
        n = h_base.sum()
        if n < 1: continue
        ax3.step(mass_cen * scale, h_base / n, where="mid",
                 color=color, lw=2.5, ls="-",
                 label=f"{label}  (×{scale:.2f})")

    ax3.axvline(16.8, color="red", lw=1.4, ls=":", label="m_X17 = 16.8 MeV")
    ax3.set_xlabel("Invariant mass  [MeV]")
    ax3.set_ylabel("Normalised counts")
    ax3.set_ylim(0, _mass_auto_ylim(
        acc.h_mass_reco[et], acc.h_mass_reco_lsbs[et], acc.h_mass_reco_all[et]))
    ax3.set_title(
        "X17 signal: invariant mass — uncorrected (dashed) vs corrected (solid)\n"
        r"Correction: rescale M by $\sqrt{1/((1-a_{e^-})(1-a_{e^+}))}$"
        "  — exact at peak for symmetric pairs; approximate in asymmetric tails")
    ax3.legend(fontsize=8, ncol=2); ax3.grid(True, alpha=0.3)
    fig3.tight_layout(); pdf.savefig(fig3); plt.close(fig3)


def plot_containment_extended(pdf, acc):
    """Extended containment QA — X17 signal only.

    Page 1 (2×3): LS+BS and all scintillators 2D containment maps + coverage/correction plots.
    Page 2 (1×4): Reco invariant mass for four energy-estimation methods.
    """
    et = 0  # X17 signal only
    ke_cen   = 0.5 * (acc.qa_ke_bins[:-1]   + acc.qa_ke_bins[1:])
    edep_cen = 0.5 * (acc.qa_edep_bins[:-1] + acc.qa_edep_bins[1:])
    _n = len(ke_cen)

    # ── Page 1: containment maps ──────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    fig.suptitle("Extended Calorimeter Containment — X17 signal only", fontsize=13)

    for row, (particle, h_ls, h_lsbs, h_all, h_bs_hit, h_reach) in enumerate([
        ("e⁻", acc.h_qa_ls_em[et], acc.h_qa_lsbs_em[et],
                acc.h_qa_all_em[et], acc.h_bs_hit_em[et], acc.h_stop_reach_em[et]),
        ("e⁺", acc.h_qa_ls_ep[et], acc.h_qa_lsbs_ep[et],
                acc.h_qa_all_ep[et], acc.h_bs_hit_ep[et], acc.h_stop_reach_ep[et]),
    ]):
        # Col 0: LS + BackScint edep, BackScint-hitting particles only
        ax = axes[row, 0]
        rs = h_lsbs.sum(axis=1, keepdims=True).clip(1)
        im = ax.pcolormesh(ke_cen, edep_cen, (h_lsbs / rs).T, cmap="viridis", vmin=0)
        ax.plot([0, 22], [0, 22], "w--", lw=1, label="perfect")
        med = _compute_medians(h_lsbs, edep_cen)
        v = ~np.isnan(med)
        if v.any():
            ax.plot(ke_cen[v], med[v], "r-o", ms=3, lw=1.5, label="median")
        plt.colorbar(im, ax=ax, label="P(LS+BS edep | E_total)")
        ax.set_xlabel("True E_total  [MeV]"); ax.set_ylabel("LS + BackScint edep  [MeV]")
        ax.set_title(f"{particle}: LS + BackScint edep\n(BackScint-hitting particles only)")
        ax.legend(fontsize=8)

        # Col 1: All scintillators edep
        ax = axes[row, 1]
        rs = h_all.sum(axis=1, keepdims=True).clip(1)
        im = ax.pcolormesh(ke_cen, edep_cen, (h_all / rs).T, cmap="viridis", vmin=0)
        ax.plot([0, 22], [0, 22], "w--", lw=1, label="perfect")
        med = _compute_medians(h_all, edep_cen)
        v = ~np.isnan(med)
        if v.any():
            ax.plot(ke_cen[v], med[v], "r-o", ms=3, lw=1.5, label="median")
        plt.colorbar(im, ax=ax, label="P(all scintillators edep | E_total)")
        ax.set_xlabel("True E_total  [MeV]"); ax.set_ylabel("All scintillators edep  [MeV]")
        ax.set_title(f"{particle}: all scintillators edep\n(PlasticScint + LS + BackScint; no MM)")
        ax.legend(fontsize=8)

    # Col 2, row 0: BackScint coverage — three categories (X17 only, e⁻+e⁺ overlaid)
    ax = axes[0, 2]
    for lbl, h_stop_ls, h_stop_bs, h_esc, h_rc, c_ls, c_bs, c_esc in [
        ("e⁻",
         acc.h_stop_stop_em[et], acc.h_stopped_in_bs_em[et],
         acc.h_escaped_ls_em[et], acc.h_stop_reach_em[et],
         "#2ca02c", "#4a90d9", "#d62728"),
        ("e⁺",
         acc.h_stop_stop_ep[et], acc.h_stopped_in_bs_ep[et],
         acc.h_escaped_ls_ep[et], acc.h_stop_reach_ep[et],
         "#76c442", "#9b59b6", "#e57373"),
    ]:
        safe = h_rc > 5
        f_ls  = np.where(safe, h_stop_ls / h_rc.clip(1), np.nan)
        f_bs  = np.where(safe, h_stop_bs / h_rc.clip(1), np.nan)
        f_esc = np.where(safe, h_esc     / h_rc.clip(1), np.nan)
        ax.plot(ke_cen[safe], f_ls[safe],  color=c_ls,  lw=2,         label=f"{lbl} stopped in LS")
        ax.plot(ke_cen[safe], f_bs[safe],  color=c_bs,  lw=2, ls="--", label=f"{lbl} stopped in BS")
        ax.plot(ke_cen[safe], f_esc[safe], color=c_esc, lw=2, ls=":",  label=f"{lbl} escaped (geo. BS miss)")
    ax.set_xlabel("True E_total  [MeV]")
    ax.set_ylabel("Fraction of LS-reaching particles")
    ax.set_title(f"BackScint coverage — 3 categories\n"
                 f"(stopped = last min-KE hit < {STOP_KE_THRESH} MeV)")
    ax.set_ylim(0, 1.05)
    ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # Col 2, row 1: Median energy estimate comparison — all methods
    ax = axes[1, 2]
    for h, lbl, col, ls in [
        (acc.h_qa_ls_em[et],   "LS only (e⁻)",        "#4a90d9", "-"),
        (acc.h_qa_lsbs_em[et], "LS+BS (e⁻, BS-hit)",  "#1a7ab5", "--"),
        (acc.h_qa_all_em[et],  "All scintillators (e⁻)",       "#0a4a85", ":"),
        (acc.h_qa_ls_ep[et],   "LS only (e⁺)",         "#e84040", "-"),
        (acc.h_qa_lsbs_ep[et], "LS+BS (e⁺, BS-hit)",  "#b82020", "--"),
        (acc.h_qa_all_ep[et],  "All scintillators (e⁺)",       "#881010", ":"),
    ]:
        med = _compute_medians(h, edep_cen)
        v = ~np.isnan(med)
        if v.any():
            ax.plot(ke_cen[v], med[v], color=col, lw=1.5, ls=ls, label=lbl)
    ax.plot([0, 22], [0, 22], "k-", lw=1.2, alpha=0.5, label="perfect (y=x)")
    ax.set_xlabel("True E_total  [MeV]")
    ax.set_ylabel("Median reconstructed energy  [MeV]")
    ax.set_title("Energy estimate comparison\n(median per method vs truth)")
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(0, 22)

    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    # ── Page 2: BackScint coverage detail — 2×2, all particles × event types ─
    # Each panel shows the three fractions (BS-hit, LS-stopped, geo-miss) vs E_total
    # for one (particle, event-type) combination.
    _LS_COL   = "#2ca02c"   # stopped in LS
    _BS_COL   = "#4a90d9"   # stopped in BS
    _MISS_COL = "#d62728"   # escaped (geometric BS miss)

    fig2, axes2 = plt.subplots(2, 2, figsize=(14, 10))
    fig2.suptitle(
        f"BackScint Coverage — Per Particle and Event Type\n"
        f"Stopped = last min-KE hit < {STOP_KE_THRESH} MeV  |  "
        f"Escaped = reached LS, not stopped in LS or BS",
        fontsize=12)

    panel_defs = [
        # ((row,col), et, particle_label, h_stop_ls, h_stop_bs, h_escaped, h_reach)
        ((0, 0), 0, "e⁻",
         acc.h_stop_stop_em[0], acc.h_stopped_in_bs_em[0],
         acc.h_escaped_ls_em[0], acc.h_stop_reach_em[0]),
        ((0, 1), 0, "e⁺",
         acc.h_stop_stop_ep[0], acc.h_stopped_in_bs_ep[0],
         acc.h_escaped_ls_ep[0], acc.h_stop_reach_ep[0]),
        ((1, 0), 1, "e⁻",
         acc.h_stop_stop_em[1], acc.h_stopped_in_bs_em[1],
         acc.h_escaped_ls_em[1], acc.h_stop_reach_em[1]),
        ((1, 1), 1, "e⁺",
         acc.h_stop_stop_ep[1], acc.h_stopped_in_bs_ep[1],
         acc.h_escaped_ls_ep[1], acc.h_stop_reach_ep[1]),
    ]
    for (r, c), event_type, particle, h_stop_ls, h_stop_bs, h_esc, h_rc in panel_defs:
        ax = axes2[r, c]
        safe = h_rc > 5
        f_ls  = np.where(safe, h_stop_ls / h_rc.clip(1), np.nan)
        f_bs  = np.where(safe, h_stop_bs / h_rc.clip(1), np.nan)
        f_esc = np.where(safe, h_esc     / h_rc.clip(1), np.nan)
        ax.fill_between(ke_cen, 0, np.where(~np.isnan(f_ls),  f_ls,  0.0), alpha=0.20, color=_LS_COL)
        ax.fill_between(ke_cen, 0, np.where(~np.isnan(f_bs),  f_bs,  0.0), alpha=0.20, color=_BS_COL)
        ax.fill_between(ke_cen, 0, np.where(~np.isnan(f_esc), f_esc, 0.0), alpha=0.20, color=_MISS_COL)
        ax.plot(ke_cen[safe], f_ls[safe],  color=_LS_COL,  lw=2.2,         label="Stopped in LS")
        ax.plot(ke_cen[safe], f_bs[safe],  color=_BS_COL,  lw=2.2, ls="--", label="Stopped in BS")
        ax.plot(ke_cen[safe], f_esc[safe], color=_MISS_COL, lw=2.2, ls=":", label="Escaped (geo. BS miss)")
        f_sum = (np.where(~np.isnan(f_ls),  f_ls,  0.0) +
                 np.where(~np.isnan(f_bs),  f_bs,  0.0) +
                 np.where(~np.isnan(f_esc), f_esc, 0.0))
        ax.plot(ke_cen[safe], f_sum[safe], color="grey", lw=0.8, ls="-",
                alpha=0.5, label="Sum (should = 1)")
        ax.set_xlabel("True E_total  [MeV]")
        ax.set_ylabel("Fraction of LS-reaching particles")
        ax.set_title(f"{ETYPE_LABEL[event_type]} — {particle}\n"
                     f"(N reach LS = {int(h_rc.sum()):,})")
        ax.set_ylim(0, 1.10)
        ax.axhline(1, color="grey", lw=0.6, ls=":")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    fig2.tight_layout(); pdf.savefig(fig2); plt.close(fig2)

    # ── Page 4: reco invariant mass — four energy methods ─────────────────
    mass_bins = acc.mass_bins
    cen       = 0.5 * (mass_bins[:-1] + mass_bins[1:])

    fig2, axes2 = plt.subplots(1, 4, figsize=(22, 5), sharey=True)
    titles = [
        "Truth mass",
        r"Reco — LS only ($E_{LS} + m_e c^2$)",
        r"Reco — LS + BackScint ($E_{LS+BS} + m_e c^2$)",
        r"Reco — all scintillators ($E_{all} + m_e c^2$)",
    ]
    hists = [
        acc.h_mass_truth[et],
        acc.h_mass_reco[et],
        acc.h_mass_reco_lsbs[et],
        acc.h_mass_reco_all[et],
    ]
    _ylim4 = _mass_auto_ylim(*hists[1:])   # from reco panels only; clips truth spike
    for ax, h, title in zip(axes2, hists, titles):
        n = h.sum()
        if n > 0:
            ax.step(cen, h / n, where="mid", color=ETYPE_COLOR[et], lw=2)
        ax.axvline(16.8, color="red", lw=1.2, ls=":", label="m_X17 = 16.8 MeV")
        ax.set_xlabel("Invariant mass  [MeV]")
        ax.set_title(title, fontsize=9)
        ax.set_ylim(0, _ylim4)
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
        ax.text(0.97, 0.96, f"N={int(n):,}", transform=ax.transAxes,
                ha="right", va="top", fontsize=8)
    axes2[0].set_ylabel("Normalised counts")
    fig2.suptitle("Reco mass — energy estimation method comparison (X17 signal)", fontsize=12)
    fig2.tight_layout(); pdf.savefig(fig2); plt.close(fig2)


def plot_trigger_vs_asymmetry(pdf, acc):
    """Trigger acceptance vs energy asymmetry — one panel per event type."""
    bins = acc.asym_bins
    cen  = 0.5 * (bins[:-1] + bins[1:])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

    trig_defs = [
        (acc.h_asym_mm_double, "MM Double — geometric ideal\n(both particles in DriftGas)", "#2ca02c", "-"),
        (acc.h_asym_single,    "Single trigger (≥1 arm)",                                   "#ff7f0e", "--"),
        (acc.h_asym_double,    "Double trigger (≥2 arms)",                                  "#1f77b4", ":"),
    ]

    for col, et in enumerate([0, 1]):
        ax = axes[col]
        tot = acc.h_asym_all[et]
        for h_dict, label, color, ls in trig_defs:
            trig = h_dict[et]
            safe = tot > 5
            frac = np.where(safe, trig / tot.clip(1), np.nan)
            ax.plot(cen[safe], frac[safe],
                    color=color, lw=2, ls=ls, marker="o", ms=3, label=label)
        ax.set_xlabel("Energy asymmetry  |KE⁻ − KE⁺| / (KE⁻ + KE⁺)")
        ax.set_title(ETYPE_LABEL[et])
        ax.set_ylim(-0.05, 1.05)
        ax.axhline(1, color="grey", lw=0.6, ls=":")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Trigger acceptance")
    fig.suptitle("Trigger acceptance vs energy asymmetry", fontsize=13)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_asymmetry(pdf, acc):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    bins = acc.asym_bins
    cen  = 0.5 * (bins[:-1] + bins[1:])

    ax = axes[0]
    for et in [0, 1]:
        h = acc.h_asym_all[et]
        norm = h.sum()
        if norm == 0: continue
        ax.step(cen, h / norm, where="mid",
                color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et], label=ETYPE_LABEL[et])
    ax.set_xlabel("Energy asymmetry  |KE⁻ − KE⁺| / (KE⁻ + KE⁺)")
    ax.set_ylabel("Normalised counts / bin")
    ax.set_title("Energy asymmetry distribution")
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1]
    for h_dict, label_sfx, ls, alpha in [
        (acc.h_asym_mm_double, "MM Double — geometric ideal", "-",   0.6),
        (acc.h_asym_single,    "single trig  (≥1 arm)",       "--",  0.8),
        (acc.h_asym_double,    "double trig  (≥2 arms)",      ":",   1.0),
    ]:
        for et in [0, 1]:
            tot  = acc.h_asym_all[et]
            trig = h_dict[et]
            safe = tot > 5
            frac = np.where(safe, trig / tot.clip(1), np.nan)
            ax.plot(cen[safe], frac[safe],
                    color=ETYPE_COLOR[et], lw=2, ls=ls, alpha=alpha, marker="o", ms=3,
                    label=f"{ETYPE_LABEL[et]}  {label_sfx}")
    ax.set_xlabel("Energy asymmetry  |KE⁻ − KE⁺| / (KE⁻ + KE⁺)")
    ax.set_ylabel("Trigger acceptance")
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.set_title("Trigger acceptance vs asymmetry\n"
                 "(drops at high asymmetry = low-KE particle not reaching trigger)")
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    fig.suptitle("Asymmetric Pair Acceptance", fontsize=13)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_pointing(pdf, acc):
    ke_bins = acc.ke_bins
    ke_cen  = 0.5 * (ke_bins[:-1] + ke_bins[1:])
    dca_cen = 0.5 * (acc.dca_bins[:-1] + acc.dca_bins[1:])  # cm

    # Distribution plot: one colour per coarse KE bin
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(ke_cen)))

    for col_idx, et in enumerate([0, 1]):
        ax = axes[col_idx]
        h = acc.h_dca[et]
        for ki in range(len(ke_cen)):
            row = h[ki]; n = row.sum()
            if n < 10: continue
            cdf = np.cumsum(row) / n
            med = dca_cen[min(np.searchsorted(cdf, 0.5), len(dca_cen)-1)]
            ax.plot(dca_cen, row / n, color=colors[ki], lw=1.5,
                    label=f"KE {ke_bins[ki]:.0f}–{ke_bins[ki+1]:.0f} MeV"
                          f"  (med={med:.2f} cm)")
        ax.set_xlabel("Track DCA to true vertex  [cm]")
        ax.set_ylabel("Normalised counts / bin")
        ax.set_title(f"MM track pointing — {ETYPE_LABEL[et]}")
        ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.3)
        ax.set_xlim(0, acc.dca_bins[-1])

    fig.suptitle("MM Track Pointing Resolution", fontsize=13)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    # Median DCA vs KE — fine binning for smooth curve
    ke_fine_cen = 0.5 * (acc.ke_fine_bins[:-1] + acc.ke_fine_bins[1:])
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    for et in [0, 1]:
        h = acc.h_dca_fine[et]; med = []
        for ki in range(len(ke_fine_cen)):
            row = h[ki]; n = row.sum()
            if n < 10: med.append(np.nan); continue
            cdf = np.cumsum(row) / n
            med.append(dca_cen[min(np.searchsorted(cdf, 0.5), len(dca_cen)-1)])
        valid = ~np.isnan(med)
        ax2.plot(ke_fine_cen[valid], np.array(med)[valid],
                 color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                 label=ETYPE_LABEL[et])
    ax2.set_xlabel("Particle KE  [MeV]")
    ax2.set_ylabel("Median DCA to true vertex  [cm]")
    ax2.set_title("Track-to-vertex DCA vs energy\n"
                  "(dominated by multiple scattering in upstream material)")
    ax2.legend(); ax2.grid(True, alpha=0.3); ax2.set_ylim(bottom=0)
    fig2.tight_layout(); pdf.savefig(fig2); plt.close(fig2)


def plot_angle_truth_vs_reco(pdf, acc):
    """2D truth vs reconstructed opening angle, one panel per event type."""
    cen = 0.5 * (acc.open2d_bins[:-1] + acc.open2d_bins[1:])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    for col, et in enumerate([0, 1]):
        ax = axes[col]
        h = acc.h_open_2d[et]
        # Normalise each truth row: P(reco | truth)
        row_sums = h.sum(axis=1, keepdims=True).clip(1)
        h_norm = h / row_sums
        im = ax.pcolormesh(cen, cen, h_norm.T, cmap="plasma", vmin=0)
        ax.plot([0, 180], [0, 180], "w--", lw=1.2, label="perfect reco")
        plt.colorbar(im, ax=ax, label="P(θ_reco | θ_truth)")
        ax.set_xlabel("Truth opening angle  [deg]")
        ax.set_title(ETYPE_LABEL[et])
        ax.set_xlim(0, 180); ax.set_ylim(0, 180)
        ax.legend(fontsize=8)

    axes[0].set_ylabel("Reconstructed opening angle  [deg]")
    fig.suptitle("Truth vs Reconstructed Opening Angle", fontsize=13)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_angle_reco(pdf, acc):
    open_cen  = 0.5 * (acc.open_bins[:-1]  + acc.open_bins[1:])
    delta_cen = 0.5 * (acc.delta_bins[:-1] + acc.delta_bins[1:])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for col_idx, et in enumerate([0, 1]):
        ax = axes[col_idx]
        h = acc.h_delta_open[et]
        col_sums = h.sum(axis=1, keepdims=True).clip(1)
        h_norm   = h / col_sums
        im = ax.pcolormesh(open_cen, delta_cen, h_norm.T,
                           cmap="plasma", vmin=0, vmax=h_norm.max())
        ax.axhline(0, color="white", lw=1, ls="--")
        plt.colorbar(im, ax=ax, label="Normalised density")
        ax.set_xlabel("Truth opening angle  [deg]")
        ax.set_ylabel("Δθ = θ_reco − θ_truth  [deg]")
        ax.set_title(f"Angle residuals — {ETYPE_LABEL[et]}")
        ax.grid(False)
    fig.suptitle("Opening Angle Reconstruction", fontsize=13)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    fig2, ax2 = plt.subplots(figsize=(9, 5))
    for et in [0, 1]:
        n   = acc.open_rms_n[et].clip(1)
        rms = np.sqrt(acc.open_rms_sum[et] / n)
        valid = acc.open_rms_n[et] > 10
        ax2.plot(open_cen[valid], rms[valid],
                 color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                 marker="o", ms=3, label=ETYPE_LABEL[et])
    ax2.set_xlabel("Truth opening angle  [deg]")
    ax2.set_ylabel("RMS(Δθ)  [deg]")
    ax2.set_title("Opening angle resolution vs truth angle\n"
                  "(RMS of θ_reco − θ_truth; dominated by MS in He-3 capsule + air)")
    ax2.legend(); ax2.grid(True, alpha=0.3); ax2.set_ylim(bottom=0)

    fig3, ax3 = plt.subplots(figsize=(8, 5))
    for et in [0, 1]:
        h = acc.h_delta_open[et].sum(axis=0); n = h.sum()
        if n == 0: continue
        rms_all = np.sqrt(acc.open_rms_sum[et].sum() / max(acc.open_rms_n[et].sum(), 1))
        ax3.step(delta_cen, h / n, where="mid",
                 color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                 label=f"{ETYPE_LABEL[et]}  RMS={rms_all:.2f}°")
    ax3.axvline(0, color="grey", lw=0.8, ls="--")
    ax3.set_xlabel("Δθ = θ_reco − θ_truth  [deg]")
    ax3.set_ylabel("Normalised counts")
    ax3.set_title("Opening angle residual (all events with MM hits)")
    ax3.legend(); ax3.grid(True, alpha=0.3)

    fig2.tight_layout(); pdf.savefig(fig2); plt.close(fig2)
    fig3.tight_layout(); pdf.savefig(fig3); plt.close(fig3)


def plot_mm_topology(pdf, acc):
    """4-page MM topology analysis: geometric classification + trigger efficiency.

    Page 1: Event topology overview (No MM / Same-arm / Diff-arm fractions)
    Page 2: Same-arm events — trigger efficiency bar + vs asymmetry
    Page 3: Diff-arm events — trigger breakdown bar + vs asymmetry
    Page 4: Combined asymmetry view — all classes normalised to all events
    """
    bins = acc.asym_bins
    cen  = 0.5 * (bins[:-1] + bins[1:])

    TOPO_COLORS = {
        "no_mm":   "#aaaaaa",
        "same":    "#ff7f0e",
        "diff":    "#2ca02c",
    }
    TRIG_COLORS = {
        "double":  "#1f77b4",
        "single":  "#ff7f0e",
        "notrig":  "#d62728",
    }

    # ── Page 1: Topology overview ─────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for col, et in enumerate([0, 1]):
        ax = axes[col]
        n = acc.n_total[et]
        if n == 0:
            continue
        n_none = n - acc.n_mm_double[et]
        n_same = acc.n_mm_same_arm[et]
        n_diff = acc.n_mm_diff_arm[et]
        labels = ["No MM hit\n(geometric or stat. miss)",
                  "Same-arm MM\n(both particles, 1 arm)",
                  "Diff-arm MM\n(particles on diff. arms)"]
        counts = [n_none, n_same, n_diff]
        colors = [TOPO_COLORS["no_mm"], TOPO_COLORS["same"], TOPO_COLORS["diff"]]
        x = np.arange(3)
        bars = ax.bar(x, [c / n for c in counts], color=colors,
                      edgecolor="k", linewidth=0.6, alpha=0.85)
        for bar, cnt in zip(bars, counts):
            frac = cnt / n
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.008,
                    f"{frac*100:.1f}%\n({cnt:,})",
                    ha="center", va="bottom", fontsize=9)
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("Fraction of all generated events")
        ax.set_title(f"{ETYPE_LABEL[et]}  (N={n:,})")
        ax.set_ylim(0, 1.15)
        ax.axhline(1, color="grey", lw=0.6, ls=":")
        ax.grid(True, axis="y", alpha=0.3)
    fig.suptitle("MM Detector Topology — Event Classification", fontsize=13)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    # ── Page 2: Same-arm trigger efficiency ──────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    for col, et in enumerate([0, 1]):
        n_same = acc.n_mm_same_arm[et]

        # Row 0: bar chart
        ax = axes[0, col]
        if n_same > 0:
            n_s   = acc.n_same_arm_single[et]
            n_nt  = acc.n_same_arm_no_trig[et]
            lbls  = ["Single trigger\n(PlScint∧LS1, ≥1 arm)", "No trigger"]
            vals  = [n_s / n_same, n_nt / n_same]
            clrs  = [TRIG_COLORS["single"], TRIG_COLORS["notrig"]]
            bars  = ax.bar(np.arange(2), vals, color=clrs,
                           edgecolor="k", linewidth=0.6, alpha=0.85)
            for bar, val, cnt in zip(bars, vals, [n_s, n_nt]):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        val + 0.01, f"{val*100:.1f}%\n({cnt:,})",
                        ha="center", va="bottom", fontsize=10)
        ax.set_xticks(np.arange(2)); ax.set_xticklabels(["Single trigger", "No trigger"])
        ax.set_ylabel("Fraction of same-arm MM events")
        ax.set_title(f"{ETYPE_LABEL[et]}\nSame-arm trigger breakdown  (N_same={n_same:,})")
        ax.set_ylim(0, 1.2); ax.axhline(1, color="grey", lw=0.6, ls=":")
        ax.grid(True, axis="y", alpha=0.3)

        # Row 1: vs asymmetry
        ax = axes[1, col]
        tot = acc.h_asym_same_arm[et]
        safe = tot > 3
        for h, lbl, clr, ls in [
            (acc.h_asym_same_arm_single[et], "Single trigger", TRIG_COLORS["single"], "-"),
            (acc.h_asym_same_arm_notrig[et], "No trigger",     TRIG_COLORS["notrig"], "--"),
        ]:
            frac = np.where(safe, h / tot.clip(1), np.nan)
            ax.plot(cen[safe], frac[safe], color=clr, lw=2, ls=ls,
                    marker="o", ms=3, label=lbl)
        ax.set_xlabel("Energy asymmetry  |KE⁻ − KE⁺| / (KE⁻ + KE⁺)")
        ax.set_ylabel("Fraction of same-arm MM events")
        ax.set_title(f"{ETYPE_LABEL[et]}: same-arm trigger efficiency vs asymmetry")
        ax.set_ylim(-0.05, 1.05)
        ax.axhline(1, color="grey", lw=0.6, ls=":")
        ax.legend(); ax.grid(True, alpha=0.3)

    fig.suptitle("Same-arm MM Events: Trigger Efficiency\n"
                 "(can only fire Singles — if not triggered, event is lost)",
                 fontsize=12)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    # ── Page 3: Diff-arm trigger breakdown ───────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    for col, et in enumerate([0, 1]):
        n_diff = acc.n_mm_diff_arm[et]

        # Row 0: bar chart
        ax = axes[0, col]
        if n_diff > 0:
            n_d  = acc.n_diff_arm_double[et]
            n_so = acc.n_diff_arm_single_only[et]
            n_nt = acc.n_diff_arm_no_trig[et]
            lbls = ["Double trigger\n(≥2 arms)", "Single-only\n(1 arm fired)", "No trigger"]
            vals = [n_d / n_diff, n_so / n_diff, n_nt / n_diff]
            clrs = [TRIG_COLORS["double"], TRIG_COLORS["single"], TRIG_COLORS["notrig"]]
            bars = ax.bar(np.arange(3), vals, color=clrs,
                          edgecolor="k", linewidth=0.6, alpha=0.85)
            for bar, val, cnt in zip(bars, vals, [n_d, n_so, n_nt]):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        val + 0.008, f"{val*100:.1f}%\n({cnt:,})",
                        ha="center", va="bottom", fontsize=9)
        ax.set_xticks(np.arange(3))
        ax.set_xticklabels(["Double trigger", "Single-only", "No trigger"])
        ax.set_ylabel("Fraction of diff-arm MM events")
        ax.set_title(f"{ETYPE_LABEL[et]}\nDiff-arm trigger breakdown  (N_diff={n_diff:,})")
        ax.set_ylim(0, 1.2); ax.axhline(1, color="grey", lw=0.6, ls=":")
        ax.grid(True, axis="y", alpha=0.3)

        # Row 1: vs asymmetry
        ax = axes[1, col]
        tot = acc.h_asym_diff_arm[et]
        safe = tot > 3
        for h, lbl, clr, ls in [
            (acc.h_asym_diff_arm_double[et],   "Double trigger",  TRIG_COLORS["double"], "-"),
            (acc.h_asym_diff_arm_singonly[et],  "Single-only",     TRIG_COLORS["single"], "--"),
            (acc.h_asym_diff_arm_notrig[et],    "No trigger",      TRIG_COLORS["notrig"], ":"),
        ]:
            frac = np.where(safe, h / tot.clip(1), np.nan)
            ax.plot(cen[safe], frac[safe], color=clr, lw=2, ls=ls,
                    marker="o", ms=3, label=lbl)
        ax.set_xlabel("Energy asymmetry  |KE⁻ − KE⁺| / (KE⁻ + KE⁺)")
        ax.set_ylabel("Fraction of diff-arm MM events")
        ax.set_title(f"{ETYPE_LABEL[et]}: diff-arm trigger breakdown vs asymmetry")
        ax.set_ylim(-0.05, 1.05)
        ax.axhline(1, color="grey", lw=0.6, ls=":")
        ax.legend(); ax.grid(True, alpha=0.3)

    fig.suptitle("Diff-arm MM Events: Trigger Breakdown\n"
                 "(ideal topology — can fire Singles or Doubles)",
                 fontsize=12)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    # ── Page 4: Combined asymmetry — all classes, norm to total events ────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
    for col, et in enumerate([0, 1]):
        ax = axes[col]
        tot = acc.h_asym_all[et]
        safe = tot > 5
        entries = [
            # (histogram, label, color, linestyle)
            (acc.h_asym_same_arm[et],          "Same-arm MM (total)",             TOPO_COLORS["same"],   "-"),
            (acc.h_asym_same_arm_single[et],   "Same-arm + Single trig",          "#d67030",             "--"),
            (acc.h_asym_same_arm_notrig[et],   "Same-arm + No trig (lost)",       "#c03000",             ":"),
            (acc.h_asym_diff_arm[et],           "Diff-arm MM (total)",             TOPO_COLORS["diff"],   "-"),
            (acc.h_asym_diff_arm_double[et],    "Diff-arm + Double trig",          TRIG_COLORS["double"], "--"),
            (acc.h_asym_diff_arm_singonly[et],  "Diff-arm + Single-only trig",     TRIG_COLORS["single"], "--"),
            (acc.h_asym_diff_arm_notrig[et],    "Diff-arm + No trig (lost)",       TRIG_COLORS["notrig"], ":"),
        ]
        for h, lbl, clr, ls in entries:
            frac = np.where(safe, h / tot.clip(1), np.nan)
            ax.plot(cen[safe], frac[safe], color=clr, lw=1.8, ls=ls,
                    marker="o", ms=2, label=lbl)
        ax.set_xlabel("Energy asymmetry  |KE⁻ − KE⁺| / (KE⁻ + KE⁺)")
        ax.set_title(ETYPE_LABEL[et])
        ax.set_ylim(-0.02, 0.85)
        ax.legend(fontsize=7.5, loc="upper right")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("Fraction of all generated events")
    fig.suptitle("MM Topology Classes vs Energy Asymmetry\n"
                 "(all classes normalised to total generated events)", fontsize=12)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def _hist1d_quantile(h, cen, q):
    """q-th quantile (0-1) of a 1D histogram; NaN if < 5 counts."""
    n = h.sum()
    if n < 5:
        return np.nan
    cdf = np.cumsum(h) / n
    return cen[min(np.searchsorted(cdf, q), len(cen) - 1)]


def plot_direction_errors(pdf, acc):
    """Per-particle direction-error decomposition (e⁻ + e⁺ combined).

    ψ = space angle between an estimator direction and the truth direction at
    the vertex (except 'drift' and 'fit error', which compare two estimators).
    This separates irreducible upstream multiple scattering from drift-volume
    scattering and from reconstruction-method error.
    """
    ke_cen  = 0.5 * (acc.ke_fine_bins[:-1] + acc.ke_fine_bins[1:])
    psi_cen = 0.5 * (acc.psi_bins[:-1] + acc.psi_bins[1:])

    def _h(name):
        return acc.h_psi[name][0] + acc.h_psi[name][1]   # X17 + IPC combined

    # ── Page 1: median ψ vs KE ────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    for name in ["first", "fit", "vline", "vcen", "nomline", "scint"]:
        med = _compute_medians(_h(name), psi_cen)
        v = ~np.isnan(med)
        if not v.any():
            continue
        ax.plot(ke_cen[v], med[v], color=PSI_COLOR[name], lw=2,
                marker="o", ms=3, label=PSI_LABEL[name])
    ax.set_xlabel("True particle KE  [MeV]")
    ax.set_ylabel("Median ψ vs truth-at-vertex  [deg]")
    ax.set_title("Direction estimators vs truth\n(lower = better measurement "
                 "of the direction at the vertex)")
    ax.set_ylim(bottom=0); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    ax = axes[1]
    for name in ["first", "last", "drift", "fit_vs_first"]:
        med = _compute_medians(_h(name), psi_cen)
        v = ~np.isnan(med)
        if not v.any():
            continue
        ax.plot(ke_cen[v], med[v], color=PSI_COLOR[name], lw=2,
                marker="o", ms=3, label=PSI_LABEL[name])
    ax.set_xlabel("True particle KE  [MeV]")
    ax.set_ylabel("Median ψ  [deg]")
    ax.set_yscale("log")
    ax.set_title("Error budget decomposition\n"
                 "upstream MS (black)  vs  drift-gas MS  vs  PCA fit error")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which="both")

    fig.suptitle("Per-Particle Direction Errors  (e⁻ + e⁺, X17 + IPC)",
                 fontsize=13)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    # ── Page 2: ψ distributions ──────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
    for ax, ke_min, title in [
        (axes[0], None, "all KE"),
        (axes[1], 8.0,  "KE > 8 MeV"),
    ]:
        for name in ["first", "fit", "vline", "vcen", "nomline", "scint"]:
            h2 = _h(name)
            rows = h2 if ke_min is None else h2[ke_cen > ke_min]
            h = rows.sum(axis=0)
            n = h.sum()
            if n < 10:
                continue
            med = _hist1d_quantile(h, psi_cen, 0.5)
            ax.step(psi_cen, h / n, where="mid", color=PSI_COLOR[name],
                    lw=1.8, label=f"{PSI_LABEL[name]}  (med={med:.1f}°)")
        ax.set_xlabel("ψ vs truth-at-vertex  [deg]")
        ax.set_title(title)
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("Normalised counts / 0.5°")
    fig.suptitle("Direction-Error Distributions  (e⁻ + e⁺, X17 + IPC)",
                 fontsize=13)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def plot_scatter_localization(pdf, acc):
    """Where along the path does the direction get destroyed?

    For each particle with a first MM hit:  t = distance from vertex to the
    hit along the truth direction, δ = perpendicular miss distance, ψ = total
    upstream direction error.  A single scatter at distance s from the vertex
    gives δ = (t − s)·tanψ, so s_eff = t − δ/tanψ localizes the scattering.
    Material positions: target wall r ≈ 15–16.4 mm, MM window at ~220 mm.
    """
    ke_cen  = 0.5 * (acc.ke_fine_bins[:-1] + acc.ke_fine_bins[1:])
    s_cen   = 0.5 * (acc.seff_bins[:-1] + acc.seff_bins[1:])
    h2 = acc.h_seff[0] + acc.h_seff[1]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    rs = h2.sum(axis=1, keepdims=True).clip(1)
    im = ax.pcolormesh(ke_cen, s_cen, (h2 / rs).T, cmap="viridis", vmin=0)
    plt.colorbar(im, ax=ax, label="P(s_eff | KE)")
    med = _compute_medians(h2, s_cen)
    v = ~np.isnan(med)
    ax.plot(ke_cen[v], med[v], "r-o", ms=3, lw=1.5, label="median")
    for y, lbl in [(16, "target wall (r≈15–16 mm)"), (220, "MM window (220 mm)")]:
        ax.axhline(y, color="white", lw=1.1, ls="--")
        ax.text(0.3, y + 8, lbl, color="white", fontsize=8)
    ax.set_xlabel("True particle KE  [MeV]")
    ax.set_ylabel("Effective scatter distance s_eff  [mm]")
    ax.set_title("Scattering localization vs KE")
    ax.legend(fontsize=8, loc="upper right")

    ax = axes[1]
    for rows, lbl, color in [
        (h2,                 "all KE",     "#1f77b4"),
        (h2[ke_cen > 8.0],   "KE > 8 MeV", "#e84040"),
    ]:
        h = rows.sum(axis=0)
        n = h.sum()
        if n < 10:
            continue
        med = _hist1d_quantile(h, s_cen, 0.5)
        ax.step(s_cen, h / n, where="mid", color=color, lw=2,
                label=f"{lbl}  (median = {med:.0f} mm)")
    ax.axvspan(15, 16.4, color="grey", alpha=0.4)
    ax.text(20, ax.get_ylim()[1] * 0.0, "", fontsize=8)
    for x, lbl in [(16, "target wall"), (220, "MM window")]:
        ax.axvline(x, color="grey", lw=1.0, ls="--")
    ax.annotate("target wall", xy=(16, 0), xytext=(30, 0.9),
                textcoords=("data", "axes fraction"), fontsize=9, color="grey")
    ax.annotate("MM window", xy=(220, 0), xytext=(230, 0.9),
                textcoords=("data", "axes fraction"), fontsize=9, color="grey")
    ax.set_xlabel("Effective scatter distance s_eff  [mm]")
    ax.set_ylabel("Normalised counts / 4 mm")
    ax.set_title("Effective scatter location\n"
                 "(single-scatter equivalent; distributed MS smears this)")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    fig.suptitle("Scattering Localization  (s_eff = t − δ/tanψ;  requires ψ > 1.5°)",
                 fontsize=13)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


# Material budget of the OLD-data geometry (what produced the analysed files):
# He-3 @ 300 bar (ρ=37.6 mg/cm³), Al wall 0.5 mm, CFRP 0.9 mm, MM face at
# 220 mm, ArIso drift gas.  (thickness [cm], X0 [cm], label)
_MS_LAYERS = [
    ("He-3 gas (300 bar, ~15 mm path)", 1.5,    94.32 / 37.6e-3,   "#4a90d9"),
    ("Al target wall 0.5 mm",           0.05,   8.897,             "#d62728"),
    ("CFRP target wall 0.9 mm",         0.09,   42.70 / 1.55,      "#ff7f0e"),
    ("Air gap ~20 cm",                  20.3,   36.62 / 1.205e-3,  "#2ca02c"),
    ("Mylar window 40 µm",              0.004,  28.54,             "#9467bd"),
    ("Kapton cathode 50 µm",            0.005,  28.58,             "#8c564b"),
    ("Cu cathode 9 µm",                 0.0009, 1.436,             "#e377c2"),
    ("Drift gas ArIso 30 mm",           3.0,    19.55 / 1.66e-3,   "#7f7f7f"),
]


def _highland(p_mev, x_over_x0):
    """Highland plane-projected θ0 [rad] for electrons (β from p)."""
    x = max(x_over_x0, 1e-12)
    e = np.sqrt(p_mev**2 + ME_MEV**2)
    beta = p_mev / e
    return (13.6 / (beta * p_mev)) * np.sqrt(x) * (1 + 0.038 * np.log(x))


def plot_ms_budget(pdf, acc):
    """Analytic Highland multiple-scattering budget vs measured upstream error.

    Layers use the OLD-data geometry (300 bar He-3, 0.9 mm CFRP, 22 cm arm
    distance) at normal incidence; oblique crossings add ~10–30%.  Measured
    points: median space angle between truth-at-vertex and the direction at
    the first DriftGas hit (upstream MS only — no reconstruction error).
    Median space angle for Gaussian MS = √(2 ln 2)·θ0 ≈ 1.177·θ0.
    """
    ke_cen  = 0.5 * (acc.ke_fine_bins[:-1] + acc.ke_fine_bins[1:])
    psi_cen = 0.5 * (acc.psi_bins[:-1] + acc.psi_bins[1:])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # ── Panel 0: x/X0 budget bars ─────────────────────────────────────────
    ax = axes[0]
    names  = [l[0] for l in _MS_LAYERS]
    xfrac  = [l[1] / l[2] for l in _MS_LAYERS]
    colors = [l[3] for l in _MS_LAYERS]
    y = np.arange(len(names))
    ax.barh(y, np.array(xfrac) * 100, color=colors, edgecolor="k", lw=0.4)
    p_ref = 10.0
    for yi, xf in zip(y, xfrac):
        th = np.degrees(_highland(p_ref, xf))
        ax.text(xf * 100 * 1.1, yi, f"θ₀={th:.2f}°", va="center", fontsize=8)
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("x / X₀  [%]")
    tot = sum(xfrac)
    ax.set_title(f"Material budget (old-data geometry)\n"
                 f"total x/X₀ = {tot*100:.2f}%   "
                 f"(θ₀ = {np.degrees(_highland(p_ref, tot)):.2f}° at p = {p_ref:.0f} MeV)")
    ax.grid(True, axis="x", alpha=0.3)

    # ── Panel 1: predicted vs measured vs KE ──────────────────────────────
    ax = axes[1]
    kes = np.linspace(0.5, 20, 100)
    ps  = np.sqrt((kes + ME_MEV)**2 - ME_MEV**2)
    for (lbl, t, x0, color) in _MS_LAYERS:
        th = np.degrees([_highland(p, t / x0) for p in ps])
        ax.plot(kes, th, color=color, lw=1.2, alpha=0.8, label=lbl)
    # Upstream total excludes the drift gas: the measured reference direction
    # is taken at the FIRST DriftGas hit, before traversing the drift volume.
    tot_up = sum(t / x0 for (lbl, t, x0, _c) in _MS_LAYERS
                 if not lbl.startswith("Drift gas"))
    th_tot = np.degrees([_highland(p, tot_up) for p in ps])
    ax.plot(kes, 1.177 * np.array(th_tot), "k-", lw=2.5,
            label="upstream total: median space angle (1.177·θ₀, no drift gas)")

    h_meas = acc.h_psi["first"][0] + acc.h_psi["first"][1]
    med = _compute_medians(h_meas, psi_cen)
    v = ~np.isnan(med)
    ax.plot(ke_cen[v], med[v], "r*", ms=9,
            label="measured: median ψ(truth, 1st-MM-hit dir)")
    ax.set_xlabel("Particle KE  [MeV]")
    ax.set_ylabel("Scattering angle  [deg]")
    ax.set_yscale("log")
    ax.set_ylim(1e-2, 1e2)
    ax.set_title("Highland prediction per layer vs measured upstream error\n"
                 "(thin: per-layer plane-projected θ₀;  thick: total median space angle)")
    ax.legend(fontsize=7, ncol=1, loc="lower left")
    ax.grid(True, alpha=0.3, which="both")

    fig.suptitle("Multiple-Scattering Budget — where the angular resolution is lost",
                 fontsize=13)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def plot_angle_method_comparison(pdf, acc):
    """Opening-angle reconstruction: method comparison.

    Page 1: Δθ residual distributions per method (X17 | IPC).
    Page 2: RMS(Δθ) vs truth opening angle per method (X17 | IPC).
    Page 3: σ68(Δθ) and median Δθ vs min(KE⁻, KE⁺)  (X17 + IPC combined).
    """
    delta_cen = 0.5 * (acc.delta_bins[:-1] + acc.delta_bins[1:])
    mo_cen    = 0.5 * (acc.mopen_bins[:-1] + acc.mopen_bins[1:])
    ke_cen    = 0.5 * (acc.ke_fine_bins[:-1] + acc.ke_fine_bins[1:])

    # ── Page 1: residual distributions ───────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
    for ax, et in zip(axes, [0, 1]):
        for mth in OPEN_METHODS:
            h = acc.h_mopen_delta[mth][et]
            n = h.sum()
            if n < 10:
                continue
            q16 = _hist1d_quantile(h, delta_cen, 0.16)
            q84 = _hist1d_quantile(h, delta_cen, 0.84)
            s68 = 0.5 * (q84 - q16)
            ax.step(delta_cen, h / n, where="mid", color=METH_COLOR[mth],
                    lw=1.8, label=f"{METH_LABEL[mth]}  (σ68={s68:.1f}°)")
        ax.axvline(0, color="grey", lw=0.8, ls="--")
        ax.set_xlabel("Δθ = θ_reco − θ_truth  [deg]")
        ax.set_title(ETYPE_LABEL[et])
        ax.legend(fontsize=7.5); ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("Normalised counts / 1°")
    fig.suptitle("Opening-Angle Residuals by Reconstruction Method", fontsize=13)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    # ── Page 2: RMS vs truth angle ───────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
    for ax, et in zip(axes, [0, 1]):
        for mth in OPEN_METHODS:
            n   = acc.mopen_rms_n[mth][et]
            rms = np.sqrt(acc.mopen_rms_sum[mth][et] / n.clip(1))
            v   = n > 20
            if not v.any():
                continue
            ax.plot(mo_cen[v], rms[v], color=METH_COLOR[mth], lw=1.8,
                    marker="o", ms=3, label=METH_LABEL[mth])
        ax.set_xlabel("Truth opening angle  [deg]")
        ax.set_title(ETYPE_LABEL[et])
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=7.5); ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("RMS(Δθ)  [deg]")
    fig.suptitle("Opening-Angle Resolution vs Truth Angle by Method", fontsize=13)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    # ── Page 3: resolution vs min pair KE ────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for mth in OPEN_METHODS:
        h2 = acc.h_mopen_minke[mth][0] + acc.h_mopen_minke[mth][1]
        s68, bias = np.full(len(ke_cen), np.nan), np.full(len(ke_cen), np.nan)
        for ki in range(len(ke_cen)):
            row = h2[ki]
            if row.sum() < 30:
                continue
            q16 = _hist1d_quantile(row, delta_cen, 0.16)
            q84 = _hist1d_quantile(row, delta_cen, 0.84)
            s68[ki]  = 0.5 * (q84 - q16)
            bias[ki] = _hist1d_quantile(row, delta_cen, 0.5)
        v = ~np.isnan(s68)
        if not v.any():
            continue
        axes[0].plot(ke_cen[v], s68[v], color=METH_COLOR[mth], lw=1.8,
                     marker="o", ms=3, label=METH_LABEL[mth])
        axes[1].plot(ke_cen[v], bias[v], color=METH_COLOR[mth], lw=1.8,
                     marker="o", ms=3, label=METH_LABEL[mth])
    axes[0].set_ylabel("σ68(Δθ)  [deg]")
    axes[0].set_title("Resolution vs min pair energy")
    axes[0].set_ylim(bottom=0)
    axes[1].set_ylabel("Median Δθ (bias)  [deg]")
    axes[1].set_title("Bias vs min pair energy")
    axes[1].axhline(0, color="grey", lw=0.8, ls="--")
    for ax in axes:
        ax.set_xlabel("min(KE⁻, KE⁺)  [MeV]")
        ax.legend(fontsize=7.5); ax.grid(True, alpha=0.3)
    fig.suptitle("Opening-Angle Performance vs Pair Energy  (X17 + IPC combined)\n"
                 "→ an energy cut directly buys angular resolution", fontsize=12)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def plot_angle_energy_explanation(pdf, acc):
    """Why is σ68(Δθ) ≈ 13° when the per-particle MS budget says θ₀ ≈ 7°?

    Three multiplicative effects:
      1. TWO particles contribute: each adds its in-decay-plane scattering
         component (one plane projection, width θ₀) → quadrature sum.
      2. For X17, KE⁻ + KE⁺ = E_trans − 2mₑ ≈ 19.6 MeV, so one particle is
         always below ~9.8 MeV and θ₀ ∝ 1/p amplifies the soft one.
      3. The quoted 6.9° is θ₀ at p = 10 MeV; the population average runs
         over the full spectrum, weighted toward lower momenta.
    """
    KE_SUM = 20.58 - 2 * ME_MEV   # X17: total pair KE [MeV]
    delta_cen = 0.5 * (acc.delta_bins[:-1] + acc.delta_bins[1:])
    ke_cen    = 0.5 * (acc.ke_fine_bins[:-1] + acc.ke_fine_bins[1:])
    tot_up = sum(t / x0 for (lbl, t, x0, _c) in _MS_LAYERS
                 if not lbl.startswith("Drift gas"))

    def _th0_deg(ke):
        p = np.sqrt((ke + ME_MEV)**2 - ME_MEV**2)
        return np.degrees(_highland(p, tot_up))

    h2 = acc.h_mopen_minke["first"][0]          # X17 signal, current method
    s68  = np.full(len(ke_cen), np.nan)
    rms  = np.full(len(ke_cen), np.nan)
    nrow = h2.sum(axis=1)
    for ki in range(len(ke_cen)):
        if nrow[ki] < 30:
            continue
        q16 = _hist1d_quantile(h2[ki], delta_cen, 0.16)
        q84 = _hist1d_quantile(h2[ki], delta_cen, 0.84)
        s68[ki] = 0.5 * (q84 - q16)
        rms[ki] = np.sqrt((h2[ki] * delta_cen**2).sum() / nrow[ki])

    # Highland prediction: quadrature of both particles' in-plane projections
    kes  = np.linspace(0.5, KE_SUM / 2, 200)
    pred = np.sqrt(_th0_deg(kes)**2 + _th0_deg(np.maximum(KE_SUM - kes, 0.1))**2)

    # Spectrum-averaged prediction, weighted by the accepted min-KE spectrum
    w  = nrow / max(nrow.sum(), 1)
    ok = (ke_cen > 0.5) & (ke_cen < KE_SUM / 2) & (w > 0)
    pred_at_bins = np.sqrt(_th0_deg(ke_cen)**2 +
                           _th0_deg(np.maximum(KE_SUM - ke_cen, 0.1))**2)
    pred_avg = np.sqrt((w[ok] * pred_at_bins[ok]**2).sum() / w[ok].sum())

    h_all = acc.h_mopen_delta["first"][0]
    q16 = _hist1d_quantile(h_all, delta_cen, 0.16)
    q84 = _hist1d_quantile(h_all, delta_cen, 0.84)
    s68_overall = 0.5 * (q84 - q16)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    v = ~np.isnan(s68)
    ax.plot(ke_cen[v], s68[v], "ko-", ms=4, lw=2, label="measured σ68(Δθ)")
    ax.plot(ke_cen[v], rms[v], "o-", ms=3, lw=1.5, color="#7f7f7f",
            label="measured RMS(Δθ)")
    ax.plot(kes, pred, "r--", lw=2,
            label=r"Highland: $\sqrt{\theta_0^2(KE_{min}) + \theta_0^2(KE_{sum}-KE_{min})}$")
    ax.axhline(s68_overall, color="k", lw=1, ls=":",
               label=f"spectrum-averaged σ68 = {s68_overall:.1f}°")
    ax.set_xlabel("min(KE⁻, KE⁺)  [MeV]")
    ax.set_ylabel("Opening-angle error  [deg]")
    ax.set_xlim(0, KE_SUM / 2 + 0.5)
    ax.set_ylim(0, 45)
    ax.set_title("X17: opening-angle residual vs pair energy\n"
                 "(current method: dir @ 1st MM hit)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.bar(ke_cen, w, width=np.diff(acc.ke_fine_bins), color="#4a90d9",
           alpha=0.6, edgecolor="none", label="accepted min-KE spectrum (X17)")
    ax.set_xlabel("min(KE⁻, KE⁺)  [MeV]")
    ax.set_ylabel("Fraction of accepted pairs")
    ax.set_xlim(0, KE_SUM / 2 + 0.5)
    txt = (
        "Per-particle budget → pair residual:\n"
        f"  θ₀(plane) at p = 10 MeV       :  {_th0_deg(10.0 - ME_MEV):.1f}°\n"
        f"  × two particles (quadrature)  :  "
        f"{np.sqrt(2) * _th0_deg(10.0 - ME_MEV):.1f}°  (at KE 9.8 + 9.8 MeV)\n"
        f"  × 1/p spectrum weighting      :  {pred_avg:.1f}°  (predicted avg)\n"
        f"  measured spectrum-avg σ68     :  {s68_overall:.1f}°\n"
        "Remainder: non-Gaussian MS tails, oblique\n"
        "wall crossings, out-of-plane 2nd-order bias"
    )
    ax.text(0.97, 0.96, txt, transform=ax.transAxes, va="top", ha="right",
            fontsize=9, fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.9))
    ax.legend(fontsize=8, loc="upper left"); ax.grid(True, alpha=0.3)
    ax.set_title("Why 13° from a 7° budget: energy spectrum weighting")

    fig.suptitle("Opening-Angle Residual vs Energy — budget reconciliation (X17)",
                 fontsize=13)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


# ── Detector-response export for the Python fast-MC ──────────────────────────

def export_response(acc, path, src_info=""):
    """Write a JSON detector-response summary for MX17_Simulation.

    Contents (all histogram tables are raw counts; consumers normalise):
      direction_error : P(ψ | KE) tables for several direction estimators,
                        ψ = space angle between estimator and truth-at-vertex.
                        Sampling a ψ from the KE row + a uniform azimuth around
                        the true direction reproduces both the opening-angle
                        resolution AND its positive bias.
      energy_response : P(edep | E_total) tables (LS-only and all-scint, e∓)
                        for particles that deposited in the LS, plus linear
                        upstream-loss fits  loss = a·E_total + b  →
                        E_est = (edep + b)/(1 − a).
      validation      : σ68/bias of Δθ vs min pair KE measured in Geant4,
                        to compare against the fast-MC after smearing.
    """
    import json
    from datetime import date

    def _i(h):
        return np.asarray(h, dtype=np.int64).tolist()

    def _fit_loss(h_ul, e_cen, l_cen):
        med   = _compute_medians(h_ul, l_cen)
        valid = (~np.isnan(med)) & (e_cen > 0.5)
        if valid.sum() < 4:
            return None
        a, b = np.polyfit(e_cen[valid], med[valid], 1)
        return [float(a), float(b)]

    e_cen = 0.5 * (acc.qa_ke_bins[:-1] + acc.qa_ke_bins[1:])
    l_cen = 0.5 * (acc.ul_loss_bins[:-1] + acc.ul_loss_bins[1:])
    delta_cen = 0.5 * (acc.delta_bins[:-1] + acc.delta_bins[1:])
    ke_cen    = 0.5 * (acc.ke_fine_bins[:-1] + acc.ke_fine_bins[1:])

    # Validation: σ68 + bias vs min pair KE (X17+IPC combined) per method
    validation = {"minke_bin_edges": acc.ke_fine_bins.tolist(), "methods": {}}
    for mth in ["first", "vline", "fit"]:
        h2 = acc.h_mopen_minke[mth][0] + acc.h_mopen_minke[mth][1]
        s68, bias = [], []
        for row in h2:
            if row.sum() < 30:
                s68.append(None); bias.append(None)
                continue
            q16 = _hist1d_quantile(row, delta_cen, 0.16)
            q84 = _hist1d_quantile(row, delta_cen, 0.84)
            s68.append(float(0.5 * (q84 - q16)))
            bias.append(float(_hist1d_quantile(row, delta_cen, 0.5)))
        validation["methods"][mth] = {"sigma68_deg": s68, "bias_deg": bias}

    out = {
        "meta": {
            "created": str(date.today()),
            "source": src_info,
            "geometry": ("OLD-data geometry: He-3 300 bar, Al 0.5 mm + CFRP 0.9 mm "
                         "target walls, MM front face 220 mm, ArIso drift gas"),
            "n_events": {ETYPE_LABEL[et]: int(acc.n_total[et]) for et in (0, 1)},
            "units": {"ke": "MeV", "psi": "deg", "edep": "MeV"},
            "notes": ("direction_error tables: rows = true particle KE bins, "
                      "cols = psi bins, values = counts (X17+IPC combined, e- and "
                      "e+ combined). energy_response tables: rows = true E_total "
                      "bins (KE + me), cols = deposited-energy bins."),
        },
        "direction_error": {
            "ke_bin_edges":  acc.ke_fine_bins.tolist(),
            "psi_bin_edges": acc.psi_bins.tolist(),
            "estimators": {
                nm: {"label": PSI_LABEL[nm],
                     "counts": _i(acc.h_psi[nm][0] + acc.h_psi[nm][1])}
                for nm in ["first", "fit", "vline", "nomline"]
            },
        },
        "energy_response": {
            "e_total_bin_edges": acc.qa_ke_bins.tolist(),
            "edep_bin_edges":    acc.qa_edep_bins.tolist(),
            "tables": {
                "ls_em":  _i(acc.h_qa_ls_em[0]  + acc.h_qa_ls_em[1]),
                "ls_ep":  _i(acc.h_qa_ls_ep[0]  + acc.h_qa_ls_ep[1]),
                "all_em": _i(acc.h_qa_all_em[0] + acc.h_qa_all_em[1]),
                "all_ep": _i(acc.h_qa_all_ep[0] + acc.h_qa_all_ep[1]),
            },
            "loss_fit_a_b": {
                "ls_em":  _fit_loss(acc.h_ul_ls_em[0]  + acc.h_ul_ls_em[1],  e_cen, l_cen),
                "ls_ep":  _fit_loss(acc.h_ul_ls_ep[0]  + acc.h_ul_ls_ep[1],  e_cen, l_cen),
                "all_em": _fit_loss(acc.h_ul_all_em[0] + acc.h_ul_all_em[1], e_cen, l_cen),
                "all_ep": _fit_loss(acc.h_ul_all_ep[0] + acc.h_ul_all_ep[1], e_cen, l_cen),
            },
        },
        "validation": validation,
    }
    with open(path, "w") as fp:
        json.dump(out, fp)
    print(f"Detector response written → {path}")


def plot_summary_page(pdf, acc):
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.axis("off")
    lines = ["MX17 Pair Simulation Analysis\n"]
    for et in [0, 1]:
        n = acc.n_total[et]
        if n == 0:
            continue
        st    = acc.n_single_trig[et]      / n * 100
        dt    = acc.n_double_trig[et]      / n * 100
        mmd   = acc.n_mm_double[et]        / n * 100
        ls1   = acc.layer_accept[et]["LiqScint_1"] / n * 100
        ls2   = acc.layer_accept[et]["LiqScint_2"] / n * 100
        n_mmd = acc.n_mm_double[et]
        # scint efficiency given MM double
        eff_s = (acc.n_mm_double_single[et] / max(n_mmd, 1)) * 100
        eff_d = (acc.n_mm_double_double[et] / max(n_mmd, 1)) * 100
        lines.append(
            f"{'─'*52}\n"
            f"{ETYPE_LABEL[et]}  ({n:,} events)\n"
            f"  MM Double (both tracks)  :  {mmd:.1f}%\n"
            f"  LiqScint_1 reached       :  {ls1:.1f}%\n"
            f"  LiqScint_2 reached       :  {ls2:.1f}%\n"
            f"  Single trigger (all evt) :  {st:.1f}%\n"
            f"  Double trigger (all evt) :  {dt:.1f}%\n"
            f"  Single | MM Double       :  {eff_s:.1f}%\n"
            f"  Double | MM Double       :  {eff_d:.1f}%\n"
        )
    ax.text(0.05, 0.95, "\n".join(lines),
            transform=ax.transAxes, va="top", ha="left",
            fontsize=11, fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="lightyellow",
                      edgecolor="grey", alpha=0.9))
    pdf.savefig(fig); plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Analyze mx17_full_sim X17+IPC pair output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("inputs", nargs="*",
                   help="ROOT files or directories (supports shell globs)")
    p.add_argument("-o", "--outfile", default="pair_analysis.pdf",
                   help="Output PDF path")
    p.add_argument("--glob", default="*.root",
                   help="File pattern when inputs are directories")
    p.add_argument("--max-files", type=int, default=None,
                   help="Cap number of files (useful for quick tests)")
    p.add_argument("--chunk-size", type=int, default=300_000,
                   help="HitTree rows per uproot chunk")
    p.add_argument("--workers", type=int,
                   default=min(4, os.cpu_count() or 1),
                   help="Parallel worker processes")
    p.add_argument("--export-response", default=None, metavar="JSON",
                   help="Write detector-response summary (direction-error + "
                        "energy tables) for the Python fast-MC to this path")
    return p.parse_args()


def collect_files(inputs, pattern):
    files = []
    for inp in inputs:
        p = Path(inp)
        if p.is_dir():
            files.extend(sorted(p.glob(pattern)))
        elif "*" in inp or "?" in inp:
            files.extend(sorted(glob.glob(inp)))
        else:
            files.append(p)
    return [str(f) for f in files if Path(f).exists()]


def main():
    args = parse_args()

    plt.rcParams.update({
        "figure.dpi": 130, "axes.titlesize": 11,
        "axes.labelsize": 10, "xtick.labelsize": 9,
        "ytick.labelsize": 9, "legend.fontsize": 9,
    })

    if not args.inputs:
        print("No input files specified. Use --help for usage.", file=sys.stderr)
        sys.exit(1)

    files = collect_files(args.inputs, args.glob)
    if not files:
        print(f"No files found matching: {args.inputs}", file=sys.stderr)
        sys.exit(1)

    if args.max_files:
        files = files[:args.max_files]

    n_workers = min(args.workers, len(files))
    print(f"Files to process : {len(files)}")
    print(f"Workers          : {n_workers}")
    print(f"Output PDF       : {args.outfile}")

    accum = Accumulator()

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {
            pool.submit(process_file, f, args.chunk_size): f
            for f in files
        }
        for future in tqdm(as_completed(futures), total=len(futures),
                           desc="Processing", unit="file"):
            fpath = futures[future]
            try:
                accum.merge(future.result())
            except Exception as e:
                print(f"\nWARNING: {Path(fpath).name}: {e}", file=sys.stderr)

    total_events = sum(accum.n_total.values())
    print(f"\nTotal events processed: {total_events:,}")
    for et in [0, 1]:
        print(f"  {ETYPE_LABEL[et]:<20}: {accum.n_total[et]:>10,}")

    print(f"\nWriting {args.outfile} ...")
    with PdfPages(args.outfile) as pdf:
        plot_summary_page(pdf, accum)
        print("  Acceptance ...")
        plot_acceptance(pdf, accum)
        print("  MM topology ...")
        plot_mm_topology(pdf, accum)
        print("  Trigger vs asymmetry ...")
        plot_trigger_vs_asymmetry(pdf, accum)
        print("  Invariant mass ...")
        plot_invariant_mass(pdf, accum)
        print("  Calorimeter QA ...")
        plot_calorimeter_qa(pdf, accum)
        print("  Containment extended ...")
        plot_containment_extended(pdf, accum)
        print("  Upstream correction ...")
        plot_upstream_correction(pdf, accum)
        print("  Asymmetry ...")
        plot_asymmetry(pdf, accum)
        print("  MM pointing ...")
        plot_pointing(pdf, accum)
        print("  Angle truth vs reco ...")
        plot_angle_truth_vs_reco(pdf, accum)
        print("  Angle reconstruction ...")
        plot_angle_reco(pdf, accum)
        print("  MS budget ...")
        plot_ms_budget(pdf, accum)
        print("  Direction errors ...")
        plot_direction_errors(pdf, accum)
        print("  Scattering localization ...")
        plot_scatter_localization(pdf, accum)
        print("  Angle method comparison ...")
        plot_angle_method_comparison(pdf, accum)
        print("  Angle vs energy explanation ...")
        plot_angle_energy_explanation(pdf, accum)

    if args.export_response:
        export_response(accum, args.export_response,
                        src_info=f"{len(files)} files: {files[0]} ...")

    print(f"Done → {args.outfile}")


if __name__ == "__main__":
    main()
