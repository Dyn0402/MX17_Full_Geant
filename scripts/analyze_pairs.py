#!/usr/bin/env python3
"""
analyze_pairs.py — MX17 Geant4 X17+IPC pair simulation analysis
=================================================================
Reads HitTree + EventTree from mx17_full_sim output and answers:

  Q1.  Acceptance     — do e± reach each detector layer?
  Q2.  Calorimetry    — does the LS stack reconstruct a useful invariant mass?
  Q3.  Asymmetry      — are low-energy particles in asymmetric pairs stopped?
  Q4.  Pointing       — do MM tracks point back to the He-3 vertex?
  Q5.  Angle reco     — how well can we reconstruct the true opening angle?

Processes multiple unmerged ROOT files without loading everything into memory.
HitTree is read in 300k-row chunks; only the columns needed are loaded.

Usage:
    python analyze_pairs.py /eos/.../pairs/x17_ipc_pairs_job*.root
    python analyze_pairs.py /eos/.../pairs/ -o analysis.pdf --max-files 20
"""

import argparse
import gc
import glob
import sys
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

    Returns a DataFrame indexed by eventID with columns:
      - em_in_<layer>, ep_in_<layer>  bool: did this particle enter the layer?
      - em_trig, ep_trig             bool: fired PlasticScint+LiqScint_1 on any arm
      - double_trig                  bool: any 2 arms each fired by any particle
      - em_edep_ls, ep_edep_ls       float: total LS edep [MeV]
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

        # LS edep
        ls_dep = p[p["detType"].isin(LS_LAYERS)].groupby("eventID")["edep"].sum() / 1e6
        out[f"{prefix}_edep_ls"] = ls_dep.reindex(out.index, fill_value=0.0)

        # First DriftGas hit per event (position + direction in world frame)
        mm = p[p["detType"] == "DriftGas"]
        if not mm.empty:
            first_mm = mm.groupby("eventID").first()[["gx","gy","gz","px","py","pz"]]
            first_mm.columns = [f"{prefix}_mm_{c}" for c in first_mm.columns]
            out = out.join(first_mm)

    # ── Double trigger: ≥2 arms where TRIGGER_LAYERS both fired by ANY particle
    both = prim[prim["detType"].isin(TRIGGER_LAYERS)]
    if not both.empty:
        arm_layers = both.groupby(["eventID", "armID"])["detType"].apply(frozenset)
        arm_trig   = arm_layers[arm_layers.apply(lambda s: TRIGGER_LAYERS.issubset(s))]
        arms_per_event = arm_trig.groupby("eventID").size()
        double_eids = arms_per_event[arms_per_event >= 2].index
        out["double_trig"] = out.index.isin(double_eids)
    else:
        out["double_trig"] = False

    return out.reset_index()


# ── Accumulator ───────────────────────────────────────────────────────────────

class Accumulator:
    """Accumulates per-file results into histograms, never storing all rows."""

    def __init__(self):
        # Q1: acceptance — counts per event type per layer
        self.n_total       = {0: 0, 1: 0}
        self.n_double_trig = {0: 0, 1: 0}
        self.n_any_mm      = {0: 0, 1: 0}
        self.layer_accept  = {0: {l: 0 for l in SCORED_LAYERS},
                              1: {l: 0 for l in SCORED_LAYERS}}

        # Q2: invariant mass reconstruction
        self.mass_bins = np.linspace(0, 25, 101)
        self.h_mass_truth = {0: np.zeros(100), 1: np.zeros(100)}
        self.h_mass_reco  = {0: np.zeros(100), 1: np.zeros(100)}

        # Q3: energy asymmetry vs acceptance
        self.asym_bins  = np.linspace(0, 1, 51)
        self.h_asym_all    = {0: np.zeros(50), 1: np.zeros(50)}
        self.h_asym_double = {0: np.zeros(50), 1: np.zeros(50)}

        # Q4: DCA of back-projected MM track to true vertex, binned by min(KE)
        self.ke_bins  = np.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 15])
        self.dca_bins = np.linspace(0, 25, 101)
        n_ke = len(self.ke_bins) - 1
        self.h_dca = {0: np.zeros((n_ke, 100)),
                      1: np.zeros((n_ke, 100))}

        # Q5: opening angle resolution
        self.open_bins  = np.linspace(60, 180, 61)
        self.delta_bins = np.linspace(-15, 15, 101)
        self.h_delta_open = {0: np.zeros((60, 100)),
                             1: np.zeros((60, 100))}
        # Also 1D RMS in truth angle bins
        self.open_rms_sum  = {0: np.zeros(60), 1: np.zeros(60)}
        self.open_rms_n    = {0: np.zeros(60), 1: np.zeros(60)}

    def update(self, merged):
        """Update accumulators from a merged (EventTree + hit summary) DataFrame."""
        for et in [0, 1]:
            m = merged[merged["event_type"] == et]
            if m.empty:
                continue
            n = len(m)
            self.n_total[et] += n

            # Q1: layer acceptance
            # After left-merge, events with no HitTree entry get NaN for all
            # hit-summary columns; fillna(False) before any boolean indexing.
            double_mask = (m.get("double_trig",
                                 pd.Series(False, index=m.index))
                           .infer_objects(copy=False).fillna(False).astype(bool))
            self.n_double_trig[et] += int(double_mask.sum())
            def _bool_col(df, col):
                s = df.get(col, pd.Series(False, index=df.index))
                return s.infer_objects(copy=False).fillna(False).astype(bool)

            mm_mask = _bool_col(m, "em_in_DriftGas") | _bool_col(m, "ep_in_DriftGas")
            self.n_any_mm[et] += int(mm_mask.sum())

            for layer in SCORED_LAYERS:
                both = _bool_col(m, f"em_in_{layer}") | _bool_col(m, f"ep_in_{layer}")
                self.layer_accept[et][layer] += int(both.sum())

            # Q2: invariant mass reconstruction
            # Truth mass
            np.add.at(self.h_mass_truth[et],
                      np.digitize(m["inv_mass"].values, self.mass_bins) - 1,
                      1)

            # Reco mass from LS edep + MM directions
            # Only for events where both particles reached LS and MM
            mm_cols_em = [f"em_mm_{c}" for c in ["px","py","pz"]]
            mm_cols_ep = [f"ep_mm_{c}" for c in ["px","py","pz"]]
            has_reco = (m.get("em_edep_ls", 0) > 0) & (m.get("ep_edep_ls", 0) > 0)
            has_mm   = all(c in m.columns for c in mm_cols_em + mm_cols_ep)
            if has_reco.any() and has_mm:
                r = m[has_reco].copy()
                E_em = r["em_edep_ls"] + ME_MEV
                E_ep = r["ep_edep_ls"] + ME_MEV
                cos_th = (r["em_mm_px"] * r["ep_mm_px"] +
                          r["em_mm_py"] * r["ep_mm_py"] +
                          r["em_mm_pz"] * r["ep_mm_pz"]).clip(-1, 1)
                mass_reco_sq = 2 * E_em * E_ep * (1 - cos_th)
                mass_reco    = np.sqrt(mass_reco_sq.clip(0))
                idx = np.clip(np.digitize(mass_reco.values, self.mass_bins) - 1,
                              0, len(self.h_mass_reco[et]) - 1)
                np.add.at(self.h_mass_reco[et], idx, 1)

            # Q3: energy asymmetry
            asym = (np.abs(m["em_ke"] - m["ep_ke"]) /
                    (m["em_ke"] + m["ep_ke"]).clip(lower=1e-9))
            idx_all = np.clip(np.digitize(asym.values, self.asym_bins) - 1,
                              0, len(self.asym_bins) - 2)
            np.add.at(self.h_asym_all[et], idx_all, 1)
            if double_mask.any():
                idx_dt = idx_all[double_mask.values]
                np.add.at(self.h_asym_double[et], idx_dt, 1)

            # Q4: DCA of back-projected MM track to true vertex
            mm_pos_em = [f"em_mm_{c}" for c in ["gx","gy","gz"]]
            mm_dir_em = [f"em_mm_{c}" for c in ["px","py","pz"]]
            mm_pos_ep = [f"ep_mm_{c}" for c in ["gx","gy","gz"]]
            mm_dir_ep = [f"ep_mm_{c}" for c in ["px","py","pz"]]

            has_mm_both = (all(c in m.columns for c in mm_pos_em + mm_dir_em + mm_pos_ep + mm_dir_ep))
            if has_mm_both:
                has_any_mm = (m[[mm_pos_em[0], mm_pos_ep[0]]].notna().any(axis=1))
                q4 = m[has_any_mm].copy()
                if not q4.empty:
                    vtx = q4[["vtx_x","vtx_y","vtx_z"]].values  # [N,3]
                    ke_min = np.minimum(q4["em_ke"].values, q4["ep_ke"].values)

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
                        P = qv[pos_cols].values       # [N,3] mm
                        d = qv[dir_cols].values        # [N,3] unit vector
                        V = qv[["vtx_x","vtx_y","vtx_z"]].values  # [N,3] mm
                        VP  = P - V                   # [N,3]
                        t   = (VP * d).sum(axis=1, keepdims=True)
                        dca = np.linalg.norm(VP - t * d, axis=1)  # [N] mm

                        ke_this  = qv[ke_col].values
                        ke_idx   = np.clip(np.digitize(ke_this, self.ke_bins) - 1,
                                           0, len(self.ke_bins) - 2)
                        dca_idx  = np.clip(np.digitize(dca, self.dca_bins) - 1,
                                           0, len(self.dca_bins) - 2)
                        for ki, di in zip(ke_idx, dca_idx):
                            self.h_dca[et][ki, di] += 1

            # Q5: opening angle reconstruction
            has_mm_q5 = has_mm_both and (all(c in m.columns
                                              for c in mm_dir_em + mm_dir_ep))
            if has_mm_q5:
                q5 = m[m[mm_dir_em[0]].notna() & m[mm_dir_ep[0]].notna()].copy()
                if not q5.empty:
                    d_em = q5[mm_dir_em].values
                    d_ep = q5[mm_dir_ep].values
                    cos_reco = (d_em * d_ep).sum(axis=1).clip(-1, 1)
                    theta_reco  = np.degrees(np.arccos(cos_reco))
                    theta_truth = q5["openingAngle"].values
                    delta       = theta_reco - theta_truth

                    truth_idx = np.clip(
                        np.digitize(theta_truth, self.open_bins) - 1,
                        0, len(self.open_bins) - 2)
                    delta_idx = np.clip(
                        np.digitize(delta, self.delta_bins) - 1,
                        0, len(self.delta_bins) - 2)

                    for ti, di in zip(truth_idx, delta_idx):
                        self.h_delta_open[et][ti, di] += 1

                    # RMS per truth-angle bin
                    np.add.at(self.open_rms_sum[et], truth_idx, delta**2)
                    np.add.at(self.open_rms_n[et],   truth_idx, 1)


# ── Per-file processing ───────────────────────────────────────────────────────

def process_file(filepath, accum, chunk_size=300_000):
    """Process one ROOT file and update the Accumulator in place."""
    with uproot.open(filepath) as f:
        if "EventTree" not in f or "HitTree" not in f:
            print(f"  WARNING: missing trees in {filepath}", file=sys.stderr)
            return 0

        evt_df = f["EventTree"].arrays(library="pd")
        if evt_df.empty:
            return 0

        # Buffer partial last event across chunks
        leftover = None
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
        n = len(merged)
        del merged
        gc.collect()
        return n


# ── Plotting ──────────────────────────────────────────────────────────────────

def _legend_patches():
    return [mpatches.Patch(color=ETYPE_COLOR[e], label=ETYPE_LABEL[e])
            for e in [0, 1]]


def plot_q1_acceptance(pdf, acc):
    """Q1: fraction of events where e± reach each layer."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: per-layer acceptance fraction (any particle hit that layer)
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
    ax.set_ylabel("Fraction of events with ≥1 hit (any particle)")
    ax.set_title("Q1: Layer acceptance (e⁻ or e⁺)")
    ax.set_ylim(0, 1.1); ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.legend(); ax.grid(True, axis="y", alpha=0.3)

    # Right: double trigger rate and any-MM rate
    ax = axes[1]
    labels = ["Any MM hit\n(track reco possible)",
              "Double trigger\n(PlScint∧LS1 on ≥2 arms)"]
    for i, et in enumerate([0, 1]):
        n = acc.n_total[et]
        if n == 0:
            continue
        vals = [acc.n_any_mm[et] / n, acc.n_double_trig[et] / n]
        ax.bar(np.array([0, 1]) + (i - 0.5) * 0.35, vals, 0.35,
               color=ETYPE_COLOR[et], alpha=0.85, label=ETYPE_LABEL[et],
               edgecolor="k", linewidth=0.4)
        for xi, vi in zip(np.array([0, 1]) + (i - 0.5) * 0.35, vals):
            ax.text(xi, vi + 0.01, f"{vi*100:.1f}%", ha="center",
                    va="bottom", fontsize=9)

    ax.set_xticks([0, 1]); ax.set_xticklabels(labels)
    ax.set_ylabel("Fraction of events")
    ax.set_title("Q1: Key trigger rates")
    ax.set_ylim(0, 1.1); ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.legend(); ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle("Q1: Detection Acceptance", fontsize=13)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_q2_invariant_mass(pdf, acc):
    """Q2: invariant mass truth distribution and calorimeter reconstruction."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    bins = acc.mass_bins
    cen  = 0.5 * (bins[:-1] + bins[1:])

    # Left: truth invariant mass
    ax = axes[0]
    for et in [0, 1]:
        n = acc.n_total[et]
        if n == 0:
            continue
        h = acc.h_mass_truth[et]
        norm = h.sum()
        ax.step(cen, h / max(norm, 1), where="mid",
                color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                label=f"{ETYPE_LABEL[et]}  (N={int(norm):,})")
    ax.axvline(16.8, color="red", lw=1.2, ls=":", label="m_X17 = 16.8 MeV")
    ax.set_xlabel("True e⁺e⁻ invariant mass  [MeV]")
    ax.set_ylabel("Normalised counts")
    ax.set_title("Q2: Truth invariant mass distribution")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.set_yscale("log"); ax.set_ylim(bottom=1e-5)

    # Right: reconstructed invariant mass from LS edep + MM directions
    ax = axes[1]
    for et in [0, 1]:
        h = acc.h_mass_reco[et]
        norm = h.sum()
        if norm == 0:
            continue
        ax.step(cen, h / norm, where="mid",
                color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                label=f"{ETYPE_LABEL[et]}  (N={int(norm):,})")
    ax.axvline(16.8, color="red", lw=1.2, ls=":", label="m_X17 = 16.8 MeV")
    ax.set_xlabel("Reconstructed invariant mass  [MeV]\n"
                  "(from LS edep + MM track directions)")
    ax.set_ylabel("Normalised counts")
    ax.set_title("Q2: Reconstructed invariant mass\n"
                 "(events where both particles reached LS)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.set_yscale("log"); ax.set_ylim(bottom=1e-5)

    fig.suptitle("Q2: Calorimeter Invariant Mass Reconstruction", fontsize=13)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_q3_asymmetry(pdf, acc):
    """Q3: energy asymmetry vs double-trigger acceptance."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    bins = acc.asym_bins
    cen  = 0.5 * (bins[:-1] + bins[1:])

    # Left: energy asymmetry distribution
    ax = axes[0]
    for et in [0, 1]:
        h = acc.h_asym_all[et]
        norm = h.sum()
        if norm == 0:
            continue
        ax.step(cen, h / norm, where="mid",
                color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                label=ETYPE_LABEL[et])
    ax.set_xlabel("Energy asymmetry  |KE⁻ − KE⁺| / (KE⁻ + KE⁺)")
    ax.set_ylabel("Normalised counts / bin")
    ax.set_title("Q3: Energy asymmetry distribution")
    ax.legend(); ax.grid(True, alpha=0.3)

    # Right: double-trigger acceptance vs asymmetry
    ax = axes[1]
    for et in [0, 1]:
        tot = acc.h_asym_all[et]
        dt  = acc.h_asym_double[et]
        safe = tot > 5
        frac = np.where(safe, dt / tot.clip(1), np.nan)
        ax.plot(cen[safe], frac[safe],
                color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                marker="o", ms=3, label=ETYPE_LABEL[et])
    ax.set_xlabel("Energy asymmetry  |KE⁻ − KE⁺| / (KE⁻ + KE⁺)")
    ax.set_ylabel("Double-trigger acceptance")
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.set_title("Q3: Double-trigger acceptance vs asymmetry\n"
                 "(low acceptance at right = one particle stopped)")
    ax.legend(); ax.grid(True, alpha=0.3)

    fig.suptitle("Q3: Asymmetric Pair Acceptance", fontsize=13)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_q4_pointing(pdf, acc):
    """Q4: DCA of back-projected MM track to true He-3 vertex."""
    ke_bins = acc.ke_bins
    ke_cen  = 0.5 * (ke_bins[:-1] + ke_bins[1:])
    dca_cen = 0.5 * (acc.dca_bins[:-1] + acc.dca_bins[1:])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(ke_cen)))

    for col_idx, et in enumerate([0, 1]):
        ax = axes[col_idx]
        h = acc.h_dca[et]  # [n_ke_bins, n_dca_bins]
        for ki in range(len(ke_cen)):
            row = h[ki]
            n = row.sum()
            if n < 10:
                continue
            # Median DCA for this KE bin
            cdf = np.cumsum(row) / n
            med_idx = np.searchsorted(cdf, 0.5)
            med_dca = dca_cen[min(med_idx, len(dca_cen)-1)]
            ax.plot(dca_cen, row / n, color=colors[ki], lw=1.5,
                    label=f"KE {ke_bins[ki]:.0f}–{ke_bins[ki+1]:.0f} MeV"
                          f"  (med={med_dca:.1f} mm)")
        ax.set_xlabel("Track DCA to true vertex  [mm]")
        ax.set_ylabel("Normalised counts / bin")
        ax.set_title(f"Q4: MM track pointing — {ETYPE_LABEL[et]}")
        ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.3)
        ax.set_xlim(0, acc.dca_bins[-1])

    # Overplot median DCA vs KE
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    for et in [0, 1]:
        h   = acc.h_dca[et]
        med = []
        for ki in range(len(ke_cen)):
            row = h[ki]
            n   = row.sum()
            if n < 10:
                med.append(np.nan)
                continue
            cdf     = np.cumsum(row) / n
            med_idx = np.searchsorted(cdf, 0.5)
            med.append(dca_cen[min(med_idx, len(dca_cen)-1)])
        valid = ~np.isnan(med)
        ax2.plot(ke_cen[valid], np.array(med)[valid],
                 color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                 marker="o", ms=5, label=ETYPE_LABEL[et])
    ax2.set_xlabel("Particle KE  [MeV]")
    ax2.set_ylabel("Median DCA to true vertex  [mm]")
    ax2.set_title("Q4: Track-to-vertex DCA vs energy\n"
                  "(dominated by multiple scattering in upstream material)")
    ax2.legend(); ax2.grid(True, alpha=0.3); ax2.set_ylim(bottom=0)

    fig.suptitle("Q4: MM Track Pointing Resolution", fontsize=13)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)
    fig2.tight_layout()
    pdf.savefig(fig2); plt.close(fig2)


def plot_q5_angle_reco(pdf, acc):
    """Q5: reconstructed opening angle vs truth."""
    open_cen  = 0.5 * (acc.open_bins[:-1]  + acc.open_bins[1:])
    delta_cen = 0.5 * (acc.delta_bins[:-1] + acc.delta_bins[1:])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: 2D heatmap of (truth angle, reco-truth angle)
    for col_idx, et in enumerate([0, 1]):
        ax = axes[col_idx]
        h = acc.h_delta_open[et]
        # Normalise each truth-angle column
        col_sums = h.sum(axis=1, keepdims=True).clip(1)
        h_norm   = h / col_sums
        im = ax.pcolormesh(open_cen, delta_cen, h_norm.T,
                           cmap="plasma", vmin=0, vmax=h_norm.max())
        ax.axhline(0, color="white", lw=1, ls="--")
        plt.colorbar(im, ax=ax, label="Normalised density")
        ax.set_xlabel("Truth opening angle  [deg]")
        ax.set_ylabel("Δθ = θ_reco − θ_truth  [deg]")
        ax.set_title(f"Q5: Angle residuals — {ETYPE_LABEL[et]}")
        ax.grid(False)

    fig.suptitle("Q5: Opening Angle Reconstruction", fontsize=13)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)

    # RMS vs truth angle
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
    ax2.set_title("Q5: Opening angle resolution vs truth angle\n"
                  "(RMS of θ_reco − θ_truth; dominated by MS in He-3 capsule + air)")
    ax2.legend(); ax2.grid(True, alpha=0.3); ax2.set_ylim(bottom=0)

    # 1D residual distribution (all truth angles combined)
    fig3, ax3 = plt.subplots(figsize=(8, 5))
    for et in [0, 1]:
        h = acc.h_delta_open[et].sum(axis=0)
        n = h.sum()
        if n == 0:
            continue
        rms_all = np.sqrt((acc.open_rms_sum[et].sum() /
                           max(acc.open_rms_n[et].sum(), 1)))
        ax3.step(delta_cen, h / n, where="mid",
                 color=ETYPE_COLOR[et], lw=2, ls=ETYPE_LS[et],
                 label=f"{ETYPE_LABEL[et]}  RMS={rms_all:.2f}°")
    ax3.axvline(0, color="grey", lw=0.8, ls="--")
    ax3.set_xlabel("Δθ = θ_reco − θ_truth  [deg]")
    ax3.set_ylabel("Normalised counts")
    ax3.set_title("Q5: Opening angle residual (all events with MM hits)")
    ax3.legend(); ax3.grid(True, alpha=0.3)

    fig2.tight_layout(); pdf.savefig(fig2); plt.close(fig2)
    fig3.tight_layout(); pdf.savefig(fig3); plt.close(fig3)


def plot_summary_page(pdf, acc):
    """Cover page with key numbers."""
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.axis("off")
    lines = ["MX17 Pair Simulation Analysis\n"]
    for et in [0, 1]:
        n = acc.n_total[et]
        if n == 0:
            continue
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

    print(f"Files to process : {len(files)}")
    print(f"Output PDF       : {args.outfile}")

    accum = Accumulator()

    for fpath in tqdm(files, desc="Processing", unit="file"):
        try:
            process_file(fpath, accum, chunk_size=args.chunk_size)
        except Exception as e:
            print(f"\nWARNING: {Path(fpath).name}: {e}", file=sys.stderr)

    # Read totals from the accumulator — reliable even if a file raised an exception.
    total_events = sum(accum.n_total.values())
    print(f"\nTotal events processed: {total_events:,}")
    for et in [0, 1]:
        print(f"  {ETYPE_LABEL[et]:<20}: {accum.n_total[et]:>10,}")

    print(f"\nWriting {args.outfile} ...")
    with PdfPages(args.outfile) as pdf:
        plot_summary_page(pdf, accum)
        print("  Q1: acceptance ...")
        plot_q1_acceptance(pdf, accum)
        print("  Q2: invariant mass ...")
        plot_q2_invariant_mass(pdf, accum)
        print("  Q3: asymmetry ...")
        plot_q3_asymmetry(pdf, accum)
        print("  Q4: MM pointing ...")
        plot_q4_pointing(pdf, accum)
        print("  Q5: angle reconstruction ...")
        plot_q5_angle_reco(pdf, accum)

    print(f"Done → {args.outfile}")


if __name__ == "__main__":
    main()
