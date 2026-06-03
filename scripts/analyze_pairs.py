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
TRIGGER_LAYERS = frozenset(["PlasticScint", "LiqScint_1"])
LS_LAYERS      = frozenset(["LiqScint_1", "LiqScint_2"])

ETYPE_LABEL = {0: "X17 signal", 1: "IPC background"}
ETYPE_COLOR = {0: "#e84040",    1: "#4a90d9"}
ETYPE_LS    = {0: "-",          1: "--"}

HIT_COLS = ["eventID", "trackID", "parentID", "armID",
            "detType", "particle", "edep",
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

def _process_hits(df):
    """
    Convert a complete-events hits DataFrame into a per-event summary.

    Returns a DataFrame with columns:
      - em_in_<layer>, ep_in_<layer>  bool: did this particle enter the layer?
      - em_trig, ep_trig             bool: PlasticScint+LiqScint_1 on any arm
      - single_trig                  bool: >=1 arm with PlScint^LS1 (any particle)
      - double_trig                  bool: >=2 arms with PlScint^LS1 (any particle)
      - em_edep_ls, ep_edep_ls       float: LS-only edep [MeV]
      - em_edep_all, ep_edep_all     float: edep in ALL scored layers [MeV]
      - em_mm_{gx,gy,gz,px,py,pz}   float: first DriftGas hit (world coords)
      - ep_mm_*                       float: same for e+
    """
    for col in ["detType", "particle"]:
        df[col] = _decode_col(df[col])

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

        # LS edep and total edep across all scored layers [MeV]
        ls_dep  = p[p["detType"].isin(LS_LAYERS)].groupby("eventID")["edep"].sum() / 1e6
        all_dep = p.groupby("eventID")["edep"].sum() / 1e6
        out[f"{prefix}_edep_ls"]  = ls_dep.reindex(out.index,  fill_value=0.0)
        out[f"{prefix}_edep_all"] = all_dep.reindex(out.index, fill_value=0.0)

        # First DriftGas hit per event (position + direction in world frame)
        mm = p[p["detType"] == "DriftGas"]
        if not mm.empty:
            first_mm = mm.groupby("eventID").first()[["gx","gy","gz","px","py","pz"]]
            first_mm.columns = [f"{prefix}_mm_{c}" for c in first_mm.columns]
            out = out.join(first_mm)

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

    return out.reset_index()


# ── Accumulator ───────────────────────────────────────────────────────────────

class Accumulator:
    """Accumulates per-file results into histograms, never storing all rows."""

    def __init__(self):
        # Acceptance counts
        self.n_total       = {0: 0, 1: 0}
        self.n_single_trig = {0: 0, 1: 0}
        self.n_double_trig = {0: 0, 1: 0}
        self.n_any_mm      = {0: 0, 1: 0}
        self.layer_accept  = {0: {l: 0 for l in SCORED_LAYERS},
                              1: {l: 0 for l in SCORED_LAYERS}}

        # Invariant mass
        self.mass_bins = np.linspace(0, 25, 101)
        self.h_mass_truth        = {0: np.zeros(100), 1: np.zeros(100)}
        self.h_mass_reco         = {0: np.zeros(100), 1: np.zeros(100)}
        self.h_mass_reco_stopped = {0: np.zeros(100), 1: np.zeros(100)}

        # Calorimeter QA: 2D histograms (true KE, edep) — 0.5 MeV bins to 22 MeV
        _n = 44
        self.qa_ke_bins   = np.linspace(0, 22, _n + 1)
        self.qa_edep_bins = np.linspace(0, 22, _n + 1)
        self.h_qa_ls_em   = {0: np.zeros((_n, _n)), 1: np.zeros((_n, _n))}
        self.h_qa_ls_ep   = {0: np.zeros((_n, _n)), 1: np.zeros((_n, _n))}
        self.h_qa_all_em  = {0: np.zeros((_n, _n)), 1: np.zeros((_n, _n))}
        self.h_qa_all_ep  = {0: np.zeros((_n, _n)), 1: np.zeros((_n, _n))}
        # Stopping fraction: n_reached_ls, n_stopped per KE bin
        self.h_stop_reach_em = {0: np.zeros(_n), 1: np.zeros(_n)}
        self.h_stop_stop_em  = {0: np.zeros(_n), 1: np.zeros(_n)}
        self.h_stop_reach_ep = {0: np.zeros(_n), 1: np.zeros(_n)}
        self.h_stop_stop_ep  = {0: np.zeros(_n), 1: np.zeros(_n)}

        # Energy asymmetry
        self.asym_bins     = np.linspace(0, 1, 51)
        self.h_asym_all    = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_single = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_double = {0: np.zeros(50), 1: np.zeros(50)}

        # DCA of back-projected MM track to true vertex, binned by KE
        self.ke_bins  = np.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 15])
        self.dca_bins = np.linspace(0, 250, 101)
        n_ke = len(self.ke_bins) - 1
        self.h_dca = {0: np.zeros((n_ke, 100)),
                      1: np.zeros((n_ke, 100))}

        # Opening angle resolution
        self.open_bins  = np.linspace(60, 180, 61)
        self.delta_bins = np.linspace(-15, 15, 101)
        self.h_delta_open  = {0: np.zeros((60, 100)), 1: np.zeros((60, 100))}
        self.open_rms_sum  = {0: np.zeros(60), 1: np.zeros(60)}
        self.open_rms_n    = {0: np.zeros(60), 1: np.zeros(60)}

    def merge(self, other):
        """Add another Accumulator into this one (used after parallel processing)."""
        for et in [0, 1]:
            self.n_total[et]       += other.n_total[et]
            self.n_single_trig[et] += other.n_single_trig[et]
            self.n_double_trig[et] += other.n_double_trig[et]
            self.n_any_mm[et]      += other.n_any_mm[et]
            for layer in SCORED_LAYERS:
                self.layer_accept[et][layer] += other.layer_accept[et][layer]
            for attr in [
                "h_mass_truth", "h_mass_reco", "h_mass_reco_stopped",
                "h_qa_ls_em", "h_qa_ls_ep", "h_qa_all_em", "h_qa_all_ep",
                "h_stop_reach_em", "h_stop_stop_em",
                "h_stop_reach_ep", "h_stop_stop_ep",
                "h_asym_all", "h_asym_single", "h_asym_double",
                "h_dca", "h_delta_open", "open_rms_sum", "open_rms_n",
            ]:
                getattr(self, attr)[et] += getattr(other, attr)[et]

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

            single_mask = _bc("single_trig")
            double_mask = _bc("double_trig")
            self.n_single_trig[et] += int(single_mask.sum())
            self.n_double_trig[et] += int(double_mask.sum())
            mm_mask = _bc("em_in_DriftGas") | _bc("ep_in_DriftGas")
            self.n_any_mm[et] += int(mm_mask.sum())

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
                E_em = sub["em_edep_ls"] + ME_MEV
                E_ep = sub["ep_edep_ls"] + ME_MEV
                cos_th = (sub["em_mm_px"]*sub["ep_mm_px"] +
                          sub["em_mm_py"]*sub["ep_mm_py"] +
                          sub["em_mm_pz"]*sub["ep_mm_pz"]).clip(-1, 1)
                return np.sqrt((2 * E_em * E_ep * (1 - cos_th)).clip(0))

            has_ls = (m.get("em_edep_ls", 0) > 0) & (m.get("ep_edep_ls", 0) > 0)
            if has_ls.any() and has_mm:
                r = m[has_ls].dropna(subset=mm_cols_em + mm_cols_ep)
                if not r.empty:
                    np.add.at(self.h_mass_reco[et],
                              np.clip(np.digitize(_reco_mass(r).values,
                                                  self.mass_bins) - 1, 0, 99), 1)

            em_stopped = (_bc("em_in_LiqScint_1") | _bc("em_in_LiqScint_2")) \
                       & ~_bc("em_in_BackScintL") & ~_bc("em_in_BackScintR")
            ep_stopped = (_bc("ep_in_LiqScint_1") | _bc("ep_in_LiqScint_2")) \
                       & ~_bc("ep_in_BackScintL") & ~_bc("ep_in_BackScintR")
            both_stopped = em_stopped & ep_stopped & has_ls
            if both_stopped.any() and has_mm:
                r2 = m[both_stopped].dropna(subset=mm_cols_em + mm_cols_ep)
                if not r2.empty:
                    np.add.at(self.h_mass_reco_stopped[et],
                              np.clip(np.digitize(_reco_mass(r2).values,
                                                  self.mass_bins) - 1, 0, 99), 1)

            # ── Calorimeter QA ────────────────────────────────────────────
            for (ke_col, edep_ls_col, edep_all_col,
                 h_ls, h_all, h_reach, h_stop, is_stopped) in [
                ("em_ke", "em_edep_ls", "em_edep_all",
                 self.h_qa_ls_em[et], self.h_qa_all_em[et],
                 self.h_stop_reach_em[et], self.h_stop_stop_em[et], em_stopped),
                ("ep_ke", "ep_edep_ls", "ep_edep_all",
                 self.h_qa_ls_ep[et], self.h_qa_all_ep[et],
                 self.h_stop_reach_ep[et], self.h_stop_stop_ep[et], ep_stopped),
            ]:
                ke      = m[ke_col].values
                ki      = np.clip(np.digitize(ke, self.qa_ke_bins) - 1,
                                  0, len(self.qa_ke_bins) - 2)
                edep_ls  = m.get(edep_ls_col,  pd.Series(0, index=m.index)).fillna(0).values
                edep_all = m.get(edep_all_col, pd.Series(0, index=m.index)).fillna(0).values
                reached  = edep_ls > 0
                if reached.any():
                    ei_ls  = np.clip(np.digitize(edep_ls[reached],  self.qa_edep_bins) - 1,
                                     0, len(self.qa_edep_bins) - 2)
                    ei_all = np.clip(np.digitize(edep_all[reached], self.qa_edep_bins) - 1,
                                     0, len(self.qa_edep_bins) - 2)
                    np.add.at(h_ls,  (ki[reached], ei_ls),  1)
                    np.add.at(h_all, (ki[reached], ei_all), 1)
                    np.add.at(h_reach, ki[reached], 1)
                np.add.at(h_stop, ki[is_stopped.values], 1)

            # ── Energy asymmetry ──────────────────────────────────────────
            asym = (np.abs(m["em_ke"] - m["ep_ke"]) /
                    (m["em_ke"] + m["ep_ke"]).clip(lower=1e-9))
            idx_all = np.clip(np.digitize(asym.values, self.asym_bins) - 1,
                              0, len(self.asym_bins) - 2)
            np.add.at(self.h_asym_all[et], idx_all, 1)
            if single_mask.any():
                np.add.at(self.h_asym_single[et], idx_all[single_mask.values], 1)
            if double_mask.any():
                np.add.at(self.h_asym_double[et], idx_all[double_mask.values], 1)

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
                        ke_this = qv[ke_col].values
                        ke_idx  = np.clip(np.digitize(ke_this, self.ke_bins) - 1,
                                          0, len(self.ke_bins) - 2)
                        dca_idx = np.clip(np.digitize(dca, self.dca_bins) - 1,
                                          0, len(self.dca_bins) - 2)
                        np.add.at(self.h_dca[et], (ke_idx, dca_idx), 1)

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
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

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

    ax = axes[1]
    trig_labels = ["Any MM hit\n(track reco possible)",
                   "Single trigger\n(PlScint∧LS1, ≥1 arm)",
                   "Double trigger\n(PlScint∧LS1, ≥2 arms)"]
    x3 = np.arange(3)
    for i, et in enumerate([0, 1]):
        n = acc.n_total[et]
        if n == 0:
            continue
        vals = [acc.n_any_mm[et] / n,
                acc.n_single_trig[et] / n,
                acc.n_double_trig[et] / n]
        ax.bar(x3 + (i - 0.5) * 0.35, vals, 0.35, color=ETYPE_COLOR[et],
               alpha=0.85, label=ETYPE_LABEL[et], edgecolor="k", linewidth=0.4)
        for xi, vi in zip(x3 + (i - 0.5) * 0.35, vals):
            ax.text(xi, vi + 0.01, f"{vi*100:.1f}%", ha="center",
                    va="bottom", fontsize=8)
    ax.set_xticks(x3); ax.set_xticklabels(trig_labels, fontsize=9)
    ax.set_ylabel("Fraction of events")
    ax.set_title("Key trigger rates")
    ax.set_ylim(0, 1.1); ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.legend(); ax.grid(True, axis="y", alpha=0.3)

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
    ax.set_yscale("log"); ax.set_ylim(bottom=1e-5)

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
    ax.set_title("Reco mass — both particles reached LS\n"
                 r"$M = \sqrt{2\,E_{e^-}^{\rm LS}\,E_{e^+}^{\rm LS}\,(1-\cos\theta_{\rm MM})}$")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.set_yscale("log"); ax.set_ylim(bottom=1e-5)

    ax = axes[2]
    for et in [0, 1]:
        h = acc.h_mass_reco_stopped[et]
        norm = h.sum()
        if norm == 0: continue
        ax.step(cen, h / norm, where="mid",
                color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                label=f"{ETYPE_LABEL[et]}  (N={int(norm):,})")
    ax.axvline(16.8, color="red", lw=1.2, ls=":", label="m_X17 = 16.8 MeV")
    ax.set_xlabel("Reco invariant mass  [MeV]")
    ax.set_title("Reco mass — both contained in LS\n(no BackScint hit)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.set_yscale("log"); ax.set_ylim(bottom=1e-5)

    fig.suptitle("Invariant Mass  (see QA page for energy containment)", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_calorimeter_qa(pdf, acc):
    """LS energy containment QA."""
    ke_cen   = 0.5 * (acc.qa_ke_bins[:-1]   + acc.qa_ke_bins[1:])
    edep_cen = 0.5 * (acc.qa_edep_bins[:-1] + acc.qa_edep_bins[1:])

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Calorimeter QA: LS energy containment  (X17 + IPC combined)", fontsize=12)

    # Top row: 2D LS edep vs true KE for e- and e+
    for col, (particle, h_ls) in enumerate([
        ("e⁻", sum(acc.h_qa_ls_em.values())),
        ("e⁺", sum(acc.h_qa_ls_ep.values())),
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
        ax.set_xlabel("True KE  [MeV]"); ax.set_ylabel("LS edep  [MeV]")
        ax.set_title(f"{particle}: LS edep vs true KE  (LiqScint_1 + LiqScint_2 only)")
        ax.legend(fontsize=8)

    # Bottom-left: total edep (all scored layers) vs true KE
    ax = axes[1, 0]
    h_tot = sum(acc.h_qa_all_em.values()) + sum(acc.h_qa_all_ep.values())
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
    ax.set_xlabel("True KE  [MeV]"); ax.set_ylabel("Total scored edep  [MeV]")
    ax.set_title("e⁻ + e⁺: total edep in all scored layers vs true KE")
    ax.legend(fontsize=8)

    # Bottom-right: LS stopping fraction vs true KE
    ax = axes[1, 1]
    for particle, h_reach, h_stop, color in [
        ("e⁻", sum(acc.h_stop_reach_em.values()), sum(acc.h_stop_stop_em.values()), "#4a90d9"),
        ("e⁺", sum(acc.h_stop_reach_ep.values()), sum(acc.h_stop_stop_ep.values()), "#e84040"),
    ]:
        safe = h_reach > 5
        frac = np.where(safe, h_stop / h_reach.clip(1), np.nan)
        ax.plot(ke_cen[safe], frac[safe], color=color, lw=2, marker="o", ms=3,
                label=f"{particle} (stopped / reached LS)")
    ax.set_xlabel("True KE  [MeV]")
    ax.set_ylabel("Fraction stopped in LS\n(reached LS, no BackScint hit)")
    ax.set_title("Containment fraction vs KE")
    ax.set_ylim(0, 1.05); ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.legend(); ax.grid(True, alpha=0.3)

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
        (acc.h_asym_single, "single trig  (≥1 arm)", "-",  0.8),
        (acc.h_asym_double, "double trig  (≥2 arms)", "--", 1.0),
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
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    fig.suptitle("Asymmetric Pair Acceptance", fontsize=13)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_pointing(pdf, acc):
    ke_bins = acc.ke_bins
    ke_cen  = 0.5 * (ke_bins[:-1] + ke_bins[1:])
    dca_cen = 0.5 * (acc.dca_bins[:-1] + acc.dca_bins[1:])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(ke_cen)))

    for col_idx, et in enumerate([0, 1]):
        ax = axes[col_idx]
        h = acc.h_dca[et]
        for ki in range(len(ke_cen)):
            row = h[ki]; n = row.sum()
            if n < 10: continue
            cdf = np.cumsum(row) / n
            med_dca = dca_cen[min(np.searchsorted(cdf, 0.5), len(dca_cen)-1)]
            ax.plot(dca_cen, row / n, color=colors[ki], lw=1.5,
                    label=f"KE {ke_bins[ki]:.0f}–{ke_bins[ki+1]:.0f} MeV"
                          f"  (med={med_dca:.1f} mm)")
        ax.set_xlabel("Track DCA to true vertex  [mm]")
        ax.set_ylabel("Normalised counts / bin")
        ax.set_title(f"MM track pointing — {ETYPE_LABEL[et]}")
        ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.3)
        ax.set_xlim(0, acc.dca_bins[-1])

    fig2, ax2 = plt.subplots(figsize=(8, 5))
    for et in [0, 1]:
        h = acc.h_dca[et]; med = []
        for ki in range(len(ke_cen)):
            row = h[ki]; n = row.sum()
            if n < 10: med.append(np.nan); continue
            cdf = np.cumsum(row) / n
            med.append(dca_cen[min(np.searchsorted(cdf, 0.5), len(dca_cen)-1)])
        valid = ~np.isnan(med)
        ax2.plot(ke_cen[valid], np.array(med)[valid],
                 color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                 marker="o", ms=5, label=ETYPE_LABEL[et])
    ax2.set_xlabel("Particle KE  [MeV]")
    ax2.set_ylabel("Median DCA to true vertex  [mm]")
    ax2.set_title("Track-to-vertex DCA vs energy\n"
                  "(dominated by multiple scattering in upstream material)")
    ax2.legend(); ax2.grid(True, alpha=0.3); ax2.set_ylim(bottom=0)

    fig.suptitle("MM Track Pointing Resolution", fontsize=13)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    fig2.tight_layout(); pdf.savefig(fig2); plt.close(fig2)


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


def plot_summary_page(pdf, acc):
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.axis("off")
    lines = ["MX17 Pair Simulation Analysis\n"]
    for et in [0, 1]:
        n = acc.n_total[et]
        if n == 0:
            continue
        st   = acc.n_single_trig[et] / n * 100
        dt   = acc.n_double_trig[et] / n * 100
        mm   = acc.n_any_mm[et]      / n * 100
        ls1  = acc.layer_accept[et]["LiqScint_1"] / n * 100
        ls2  = acc.layer_accept[et]["LiqScint_2"] / n * 100
        lines.append(
            f"{'─'*48}\n"
            f"{ETYPE_LABEL[et]}  ({n:,} events)\n"
            f"  Any MM hit           :  {mm:.1f}%\n"
            f"  LiqScint_1 reached   :  {ls1:.1f}%\n"
            f"  LiqScint_2 reached   :  {ls2:.1f}%\n"
            f"  Single trigger       :  {st:.1f}%\n"
            f"  Double trigger       :  {dt:.1f}%\n"
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
        print("  Invariant mass ...")
        plot_invariant_mass(pdf, accum)
        print("  Calorimeter QA ...")
        plot_calorimeter_qa(pdf, accum)
        print("  Asymmetry ...")
        plot_asymmetry(pdf, accum)
        print("  MM pointing ...")
        plot_pointing(pdf, accum)
        print("  Angle reconstruction ...")
        plot_angle_reco(pdf, accum)

    print(f"Done → {args.outfile}")


if __name__ == "__main__":
    main()
