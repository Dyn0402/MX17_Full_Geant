#!/usr/bin/env python3
"""
analyze_full_experiment.py
Analysis of full-experiment electron/positron transmission and calorimetry.

Color scheme: same detector layer → same color (grouped by subsystem).
Line style  : solid = electron, dashed = positron.

Each transmission/edep plot is produced twice:
  - Full version  : all 10 layers
  - Simple version: MM drift gas, Plastic scint., Liq. scint. 1, Liq. scint. 2

Usage
    python3 scripts/analyze_full_experiment.py \\
        --indir /eos/user/d/dneff/mx17_geant_sim_results/full \\
        --gas ArIso \\
        --particles electron positron \\
        --outfile full_experiment_analysis.pdf
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from pathlib import Path

# Use 'spawn' start method for all worker pools.
#
# On Linux the default is 'fork', which copies the parent's memory including
# any threading locks held by matplotlib, uproot, or Python internals at the
# moment of fork().  Worker processes then inherit those locks in a locked
# state and deadlock on the IPC queue semaphore.  'spawn' starts each worker
# with a fresh Python interpreter, eliminating the inherited-lock problem at
# the cost of ~1 s per-worker startup (negligible for our batch sizes).
_MP_CTX = get_context("spawn")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import uproot
import pandas as pd

# ── Color scheme (grouped by detector subsystem) ─────────────────────────────
# He-3 target:   blue
# MM elements:   green family  (dark → light)
# PCB:           orange
# Plastic scint: red
# LS elements:   purple family (dark → light)

LAYERS = [
    # (edep_branch, primIn_branch, label, color)
    ("edepHe3Gas",      "primInHe3Gas",    "He-3 gas",        "#1f77b4"),
    ("edepDrift",       "primInDrift",     "MM drift gas",    "#2ca02c"),
    ("edepAmp",         "primInAmp",       "MM amp gas",      "#74c476"),
    ("edepResistPaste", None,              "Resistive paste", "#a1d99b"),
    ("edepPCB",         "primInPCB",       "PCB stack",       "#ff7f0e"),
    ("edepScintWall",   "primInScintWall", "Plastic scint.",  "#d62728"),
    ("edepLS1",         "primInLS1",       "Liq. scint. 1",   "#6a1d8a"),
    ("edepLS2",         "primInLS2",       "Liq. scint. 2",   "#9b59b6"),
    ("edepLS3",         "primInLS3",       "Liq. scint. 3",   "#bb8fce"),
    ("edepLS4",         "primInLS4",       "Liq. scint. 4",   "#d7bde2"),
]

# Fully granular layer list — every individually scored material.
# Used in detailed edep plots and for material budget analysis.
# Aggregation into groups is done in Python via COARSE_GROUPS.
DETAILED_LAYERS = [
    # He-3 target
    ("edepHe3Gas",       None,              "He-3 gas",          "#1f77b4"),
    # MM entrance / dead layers
    ("edepMylar",        None,              "MM Mylar window",   "#aec7e8"),
    ("edepCathode",      None,              "MM cathode (Kap+Cu)","#9edae5"),
    ("edepDrift",        "primInDrift",     "MM drift gas",      "#2ca02c"),
    ("edepMicromesh",    None,              "MM micromesh",      "#969696"),
    ("edepAmp",          "primInAmp",       "MM amp gas",        "#74c476"),
    ("edepResistPaste",  None,              "Resistive paste",   "#a1d99b"),
    # PCB individual
    ("edepPCBKapton",    None,              "PCB Kapton",        "#ffbb78"),
    ("edepPCBCu",        "primInPCB",       "PCB copper (×4)",   "#d55e00"),
    ("edepPCBFR4",       None,              "PCB FR4 (×4)",      "#ff7f0e"),
    ("edepPCBRohacell",  None,              "PCB Rohacell",      "#ffd700"),
    ("edepPCBAlFoil",    None,              "PCB Al foil",       "#c7c7c7"),
    # Scint wall
    ("edepScintTape",    None,              "Scint. tape (×2)",  "#f7b6d2"),
    ("edepScintWall",    "primInScintWall", "Plastic scint.",    "#d62728"),
    ("edepScintAlFoil",  None,              "Scint. Al foil",    "#e8a09a"),
    # LS layers
    ("edepLS1",          "primInLS1",       "Liq. scint. 1",     "#6a1d8a"),
    ("edepLS2",          "primInLS2",       "Liq. scint. 2",     "#9b59b6"),
    ("edepLS3",          "primInLS3",       "Liq. scint. 3",     "#bb8fce"),
    ("edepLS4",          "primInLS4",       "Liq. scint. 4",     "#d7bde2"),
    ("edepLSCFRP",       None,              "LS CFRP walls",     "#555555"),
]

SIMPLIFIED_LAYERS = [
    ("edepDrift",     "primInDrift",     "MM drift gas",   "#2ca02c"),
    ("edepScintWall", "primInScintWall", "Plastic scint.", "#d62728"),
    ("edepLS1",       "primInLS1",       "Liq. scint. 1",  "#6a1d8a"),
    ("edepLS4",       "primInLS4",       "Liq. scint. 4",  "#d7bde2"),
]

# Coarse grouping: aggregate detector subsystems
# (label, color, [edep branches to sum], trans_col for entry transmission)
COARSE_GROUPS = [
    ("He-3 target",        "#1f77b4",
     ["edepHe3Gas"],
     "primInHe3Gas"),
    ("Micromegas (all)",   "#2ca02c",
     ["edepDrift", "edepAmp", "edepResistPaste"],
     "primInDrift"),
    ("PCB stack",          "#ff7f0e",
     ["edepPCB"],
     "primInPCB"),
    ("Plastic scint.",     "#d62728",
     ["edepScintWall"],
     "primInScintWall"),
    ("Liquid scint. (all)","#6a1d8a",
     ["edepLS1", "edepLS2", "edepLS3", "edepLS4"],
     "primInLS1"),
]

# Cumulative path from gun (He-3 centre, z = 0) to entry face of each
# scored volume [mm].  Used for survival-vs-path plots.
# He-3 gas radius 25 mm + exit wall (Al 0.5 + CFRP 0.9) + 200 mm air
# + MM dead layers (~0.1 mm) → drift gas entry ≈ 226.5 mm.
SURVIVAL_COLS = [
    # (trans_col,        z_entry_mm, label)
    ("primInHe3Gas",       25.0,   "He-3 exit"),
    ("primInDrift",       226.5,   "MM drift"),
    ("primInAmp",         256.5,   "MM amp"),
    ("primInPCB",         256.8,   "PCB"),
    ("primInScintWall",   282.6,   "Plastic scint."),
    ("primInLS1",         307.3,   "LS1"),
    ("primInLS2",         323.8,   "LS2"),
    ("primInLS3",         340.3,   "LS3"),
    ("primInLS4",         356.8,   "LS4"),
]
LS4_EXIT_MM = 371.8   # end of LS4 (LS_CFRP_5 skipped — not scored)

SURVIVAL_ENERGIES_MEV = [1, 3, 6, 10, 16]   # integer MeV — exact points guaranteed by submit script

# Schematic stack for the diagram strip.
# (z_start_mm, z_end_mm, face_color, text_color, label)
# Air gaps are omitted — they take the figure background color.
DETECTOR_STACK = [
    ( 0.0,   25.0,  "#1f77b4", "white",  "He-3 gas"),   # He-3 gas (half-diameter)
    (25.0,   25.5,  "#aec7e8", "black",  ""),            # Al capsule wall
    (25.5,   26.4,  "#555555", "white",  ""),            # CFRP capsule wall
    # 26.4 – 226.4  Air gap 1 (background)
    (226.4,  226.5, "#c7e9c0", "black",  ""),            # MM dead layers (Mylar+cathode)
    (226.5,  256.5, "#2ca02c", "white",  "Drift"),       # Drift gas  30 mm
    (256.5,  256.55,"#808080", "white",  ""),            # Micromesh
    (256.55, 257.0, "#74c476", "black",  ""),            # Amp gas + resistive paste
    (257.0,  262.6, "#ff7f0e", "black",  "PCB"),         # PCB stack  5.6 mm
    # 262.6 – 282.6  Air gap 2 (background)
    (282.6,  285.8, "#d62728", "white",  "Scint."),      # Plastic scint. wall  3.2 mm
    # 285.8 – 305.8  Air gap 3 (background)
    (305.8,  307.3, "#555555", "white",  ""),            # CFRP  1.5 mm
    (307.3,  322.3, "#6a1d8a", "white",  "LS 1"),        # Liq. scint. 1  15 mm
    (322.3,  323.8, "#555555", "white",  ""),            # CFRP
    (323.8,  338.8, "#9b59b6", "white",  "LS 2"),        # Liq. scint. 2
    (338.8,  340.3, "#555555", "white",  ""),            # CFRP
    (340.3,  355.3, "#bb8fce", "black",  "LS 3"),        # Liq. scint. 3
    (355.3,  356.8, "#555555", "white",  ""),            # CFRP
    (356.8,  371.8, "#d7bde2", "black",  "LS 4"),        # Liq. scint. 4
]

# Line style per particle
PARTICLE_STYLES = {
    "electron": {"ls": "-",  "lw": 1.8, "label": "e⁻"},
    "positron": {"ls": "--", "lw": 1.8, "label": "e⁺"},
}

LS_EDEP_BRANCHES  = ["edepLS1", "edepLS2", "edepLS3", "edepLS4"]
SNAPSHOT_ENERGIES = [5, 8, 11, 13, 16]   # integer MeV — exact points guaranteed by submit script
MAX_EVENTS_ANG    = 3000


# ── File I/O ──────────────────────────────────────────────────────────────────

def parse_energy(tag: str):
    m = re.search(r"_(\d+(?:p\d+)?)MeV", tag)
    return float(m.group(1).replace("p", ".")) if m else None


def group_files(indir: str, gas: str, particle: str,
                prefix: str = "full") -> dict:
    """energy_MeV → [root file paths].  prefix selects 'full' or 'sr90' runs."""
    groups = defaultdict(list)
    for f in sorted(Path(indir).glob(f"{prefix}_{gas}_{particle}_*.root")):
        base   = re.sub(r"_t\d+$", "", f.stem)
        energy = parse_energy(base)
        if energy is not None:
            groups[energy].append(str(f))
    if not groups:
        print(f"  No files found for {prefix}/{gas}/{particle} in {indir}")
    return dict(sorted(groups.items()))


def read_event_tree(files: list) -> pd.DataFrame:
    dfs = []
    for fpath in files:
        try:
            with uproot.open(fpath) as f:
                if "EventTree" in f:
                    dfs.append(f["EventTree"].arrays(library="pd"))
        except Exception as e:
            print(f"  Warning: {fpath}: {e}")
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def read_cluster_tree(files: list, max_events: int = MAX_EVENTS_ANG) -> pd.DataFrame:
    """
    Read DriftGas primary-track clusters for angular analysis.

    Key performance constraints:
      entry_stop  — read at most this many rows from the tree.  The ClusterTree
                    contains entries from many volumes; DriftGas entries are
                    roughly 1 in 5.  Reading max_events×1500 rows gives a large
                    safety margin while avoiding loading the entire file (which
                    can be 15 M+ rows and take hours on EOS).
      Vectorised  — all filtering is done with numpy masks, no Python loops over
                    individual rows.
    """
    # Upper bound on rows to read: ~1 500 total entries per event (generous)
    entry_stop = max_events * 1500
    collected  = []

    for fpath in files:
        n_so_far = (int(np.unique(np.vstack(collected)[:, 0]).shape[0])
                    if collected else 0)
        if n_so_far >= max_events:
            break
        try:
            with uproot.open(fpath) as f:
                if "ClusterTree" not in f:
                    continue
                t    = f["ClusterTree"]
                stop = min(t.num_entries, entry_stop)

                a = t.arrays(["eventID", "trackID", "x", "y", "z", "volume"],
                              entry_stop=stop, library="np")

                # Decode volume column (fixed-length char array → str)
                vols = a["volume"]
                if vols.dtype.kind == "S":
                    # Strip null bytes vectorised
                    vol_str = np.frompyfunc(
                        lambda v: v.split(b"\x00")[0].decode("ascii"), 1, 1
                    )(vols).astype(str)
                else:
                    vol_str = vols.astype(str)

                # Filter: DriftGas only, primary particle (trackID == 1)
                mask = (vol_str == "DriftGas") & (a["trackID"] == 1)
                if not mask.any():
                    continue

                ev_f = a["eventID"][mask]
                x_f  = a["x"][mask]
                y_f  = a["y"][mask]
                z_f  = a["z"][mask]

                # Keep only the first max_events unique events
                unique_ev = np.unique(ev_f)
                n_needed  = max_events - n_so_far
                if len(unique_ev) > n_needed:
                    unique_ev = unique_ev[:n_needed]

                keep = np.isin(ev_f, unique_ev)
                if keep.any():
                    collected.append(np.column_stack(
                        [ev_f[keep], x_f[keep], y_f[keep], z_f[keep]]))

        except Exception as e:
            print(f"  Warning (clusters): {fpath}: {e}")

    if not collected:
        return pd.DataFrame(columns=["eventID", "x", "y", "z"])

    arr = np.vstack(collected)
    return pd.DataFrame(arr, columns=["eventID", "x", "y", "z"])


# ── Summary computation ───────────────────────────────────────────────────────

# ── Module-level worker functions (must be picklable for multiprocessing) ─────

def _summarise_one(args):
    """Process a single energy point for compute_summary."""
    energy, files = args
    df = read_event_tree(files)
    if df.empty:
        return None
    row = {"energy_MeV": energy, "n_events": len(df)}
    for (edep_br, primin_br, _label, _) in LAYERS:
        if edep_br in df.columns:
            row[f"edep_{edep_br}"]     = df[edep_br].mean() / 1e6
            row[f"edep_std_{edep_br}"] = df[edep_br].std()  / 1e6
        if primin_br and primin_br in df.columns:
            row[f"trans_{primin_br}"] = df[primin_br].mean()
    ls_cols = [b for b in LS_EDEP_BRANCHES if b in df.columns]
    if ls_cols:
        ls_tot = df[ls_cols].sum(axis=1) / 1e6
        row["ls_total_mean"]   = ls_tot.mean()
        row["ls_total_std"]    = ls_tot.std()
        row["ls_total_median"] = ls_tot.median()
        row["frac_any_ls4"]    = (df["primInLS4"].mean()
                                  if "primInLS4" in df.columns else np.nan)
        if "primInLS1" in df.columns:
            ls_in = ls_tot[df["primInLS1"].astype(bool)]
            if len(ls_in) > 10:
                row["ls_res_mean"] = ls_in.mean()
                row["ls_res_std"]  = ls_in.std()
    if "edepLSCFRP" in df.columns:
        row["edep_edepLSCFRP"] = df["edepLSCFRP"].mean() / 1e6
    if "primInLSCFRP5" in df.columns:
        row["trans_primInLSCFRP5"] = df["primInLSCFRP5"].mean()
    return row


def _trigger_one(args):
    """Process a single energy point for compute_trigger_fractions."""
    energy, files, thresh_scint_eV, thresh_ls1_eV = args
    total = fired = 0
    for fpath in files:
        try:
            with uproot.open(fpath) as f:
                if "EventTree" not in f:
                    continue
                keys = list(f["EventTree"].keys())
                need = [b for b in ("edepScintWall", "edepLS1") if b in keys]
                if len(need) < 2:
                    continue
                arr  = f["EventTree"].arrays(need, library="np")
                mask = ((arr["edepScintWall"] > thresh_scint_eV) &
                        (arr["edepLS1"]       > thresh_ls1_eV))
                fired += int(mask.sum())
                total += len(arr["edepScintWall"])
        except Exception:
            pass
    return energy, (fired / total if total > 0 else np.nan)


def _angular_one(args):
    """Load cluster data and compute angles for a single energy point."""
    energy, files = args
    cdf    = read_cluster_tree(files, MAX_EVENTS_ANG)
    angles = compute_angles(cdf)
    return energy, angles   # return full array; caller computes RMS or plots


def compute_summary(groups: dict, n_workers: int = None) -> pd.DataFrame:
    items = list(groups.items())
    n     = len(items)
    w     = n_workers or os.cpu_count()
    print(f"  {n} energy points — {min(w, n)} parallel workers")
    records = []
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=_MP_CTX) as ex:
        futures = {ex.submit(_summarise_one, item): item[0] for item in items}
        done = 0
        for fut in as_completed(futures):
            done += 1
            try:
                row = fut.result()
            except Exception as e:
                print(f"  Warning [{done}/{n}]: {e}")
                row = None
            if row is not None:
                records.append(row)
            if done % max(1, n // 10) == 0 or done == n:
                print(f"  ... {done}/{n} done")
    return pd.DataFrame(records).sort_values("energy_MeV").reset_index(drop=True)


def compute_angles(cdf: pd.DataFrame) -> np.ndarray:
    angles = []
    for _, grp in cdf.groupby("eventID"):
        if len(grp) < 5:
            continue
        z, x, y = grp["z"].values, grp["x"].values, grp["y"].values
        if z.max() - z.min() < 0.5:
            continue
        ax_ = np.polyfit(z, x, 1)[0]
        ay_ = np.polyfit(z, y, 1)[0]
        angles.append(np.degrees(np.arctan(np.sqrt(ax_**2 + ay_**2))))
    return np.array(angles)


def find_snapshots(groups: dict, targets: list) -> dict:
    avail = np.array(sorted(groups.keys()))
    return {t: (avail[np.argmin(np.abs(avail - t))],
                groups[avail[np.argmin(np.abs(avail - t))]]) for t in targets}


# ── Legend helpers ────────────────────────────────────────────────────────────

def _layer_legend(ax, layer_list, loc="lower right", ncol=2):
    """Color-patch legend for layers."""
    handles = [mpatches.Patch(color=color, label=label)
               for (_, _, label, color) in layer_list]
    return ax.legend(handles=handles, loc=loc, fontsize=7,
                     ncol=ncol, title="Layer", title_fontsize=7)


def _particle_legend(ax, summaries, loc="upper right"):
    """Line-style legend for particles."""
    handles = [mlines.Line2D([], [], color="k",
                             ls=PARTICLE_STYLES[p]["ls"],
                             lw=PARTICLE_STYLES[p]["lw"],
                             label=PARTICLE_STYLES[p]["label"])
               for p in summaries if p in PARTICLE_STYLES]
    if handles:
        return ax.legend(handles=handles, loc=loc, fontsize=9,
                         title="Particle", title_fontsize=8)
    return None


def _add_both_legends(ax, layer_list, summaries,
                      layer_loc="lower right", part_loc="upper right"):
    leg1 = _layer_legend(ax, layer_list, loc=layer_loc)
    ax.add_artist(leg1)
    _particle_legend(ax, summaries, loc=part_loc)


E_LO, E_HI = 4.0, 16.5  # physical e0 range [MeV]

def _add_energy_lines(ax):
    """Mark the physical electron energy range with vertical lines."""
    ax.axvline(E_LO, color="grey", lw=1.0, ls="--", zorder=0)
    ax.axvline(E_HI, color="grey", lw=1.0, ls="--", zorder=0)
    # Shade the region between the two lines
    ylim = ax.get_ylim()
    ax.axvspan(E_LO, E_HI, alpha=0.06, color="grey", zorder=0)
    ax.set_ylim(ylim)  # restore in case axvspan shifted it

def _set_xaxis(ax, xlabel="Electron energy (MeV)"):
    ax.set_xlabel(xlabel, fontsize=11)
    _add_energy_lines(ax)


# ── Plot functions ────────────────────────────────────────────────────────────

def _active_trans_cols(summaries: dict) -> set:
    """
    Return the set of transmission column names that have at least one
    non-zero value across all particles and energies.  Columns that are
    uniformly zero indicate a volume that does not exist in the geometry
    (e.g. He3Gas in the sr90 calibration mode) and should be hidden.
    """
    active = set()
    for summary in summaries.values():
        for col in summary.columns:
            if col.startswith("trans_") and summary[col].max() > 1e-6:
                active.add(col)
    return active


def _has_target(summaries: dict) -> bool:
    """True if the He-3 target volume was present in the simulation."""
    return "trans_primInHe3Gas" in _active_trans_cols(summaries)


def _plot_transmission_on(ax, summaries: dict, layer_list: list):
    active = _active_trans_cols(summaries)
    for (_, primin_br, label, color) in layer_list:
        if not primin_br:
            continue
        col = f"trans_{primin_br}"
        if col not in active:
            continue   # volume absent in this geometry (e.g. He-3 in sr90 mode)
        for particle, summary in summaries.items():
            if col not in summary.columns:
                continue
            pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6})
            ax.plot(summary["energy_MeV"].values, summary[col].values,
                    color=color, ls=pst["ls"], lw=pst["lw"])
    ax.set_ylim(-0.02, 1.05)
    ax.axhline(1.0, color="grey", lw=0.6, ls=":")
    ax.set_ylabel("Transmission fraction", fontsize=11)
    ax.grid(True, alpha=0.3)


def plot_transmission(pdf: PdfPages, summaries: dict):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    _plot_transmission_on(ax, summaries, LAYERS)
    _set_xaxis(ax)
    ax.set_title("Fraction reaching each layer — all layers", fontsize=12)
    _add_both_legends(ax, LAYERS, summaries)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_transmission_simple(pdf: PdfPages, summaries: dict):
    fig, ax = plt.subplots(figsize=(8, 5))
    _plot_transmission_on(ax, summaries, SIMPLIFIED_LAYERS)
    _set_xaxis(ax)
    ax.set_title("Fraction reaching each layer — key layers", fontsize=12)
    _add_both_legends(ax, SIMPLIFIED_LAYERS, summaries,
                      layer_loc="lower right", part_loc="center right")
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def _plot_edep_on(ax, summaries: dict, layer_list: list):
    for (edep_br, _, label, color) in layer_list:
        col = f"edep_{edep_br}"
        for particle, summary in summaries.items():
            if col not in summary.columns:
                continue
            vals = summary[col].values
            if np.nanmax(vals) < 1e-8:
                continue
            pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6})
            ax.plot(summary["energy_MeV"].values, vals,
                    color=color, ls=pst["ls"], lw=pst["lw"])
    ax.set_yscale("log")
    ax.set_ylabel("Mean energy deposition (MeV / primary)", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)


def plot_edep(pdf: PdfPages, summaries: dict):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    _plot_edep_on(ax, summaries, LAYERS)
    _set_xaxis(ax)
    ax.set_title("Mean energy deposited per layer (all events) — all layers",
                 fontsize=12)
    _add_both_legends(ax, LAYERS, summaries, layer_loc="upper left")
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_edep_simple(pdf: PdfPages, summaries: dict):
    fig, ax = plt.subplots(figsize=(8, 5))
    _plot_edep_on(ax, summaries, SIMPLIFIED_LAYERS)
    _set_xaxis(ax)
    ax.set_title("Mean energy deposited per layer (all events) — key layers",
                 fontsize=12)
    _add_both_legends(ax, SIMPLIFIED_LAYERS, summaries, layer_loc="upper left")
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_calorimetry(pdf: PdfPages, summaries: dict):
    if not any("ls_total_mean" in s.columns for s in summaries.values()):
        return
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: LS edep vs initial energy
    ax = axes[0]
    e_max = 0.0
    for particle, summary in summaries.items():
        if "ls_total_mean" not in summary.columns:
            continue
        pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6, "label": particle})
        en   = summary["energy_MeV"].values
        mean = summary["ls_total_mean"].values
        std  = summary["ls_total_std"].values
        ax.errorbar(en, mean, yerr=std, fmt="o", ms=2.5, lw=1,
                    ls=pst["ls"], color="#9467bd" if particle == "electron" else "#e377c2",
                    ecolor="lightgrey", capsize=1.5, label=pst["label"])
        e_max = max(e_max, en.max())
    ax.plot([0, e_max], [0, e_max], "k--", lw=0.8, label="E_dep = E₀", zorder=0)
    _add_energy_lines(ax)
    ax.set_xlabel("Initial energy (MeV)", fontsize=10)
    ax.set_ylabel("Total LS edep (MeV)", fontsize=10)
    ax.set_title("LS calorimeter response", fontsize=11)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Panel B: linearity
    ax = axes[1]
    for particle, summary in summaries.items():
        if "ls_total_mean" not in summary.columns:
            continue
        pst  = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6, "label": particle})
        en   = summary["energy_MeV"].values
        lin  = summary["ls_total_mean"].values / en
        color = "#9467bd" if particle == "electron" else "#e377c2"
        ax.plot(en, lin, ls=pst["ls"], lw=pst["lw"], color=color,
                label=pst["label"])
    ax.axhline(1.0, color="k", lw=0.8, ls=":")
    ax.set_ylim(0, 1.15)
    _add_energy_lines(ax)
    ax.set_xlabel("Initial energy (MeV)", fontsize=10)
    ax.set_ylabel("E_dep(LS) / E₀", fontsize=10)
    ax.set_title("Calorimeter linearity", fontsize=11)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Panel C: energy resolution
    ax = axes[2]
    for particle, summary in summaries.items():
        if "ls_res_std" not in summary.columns:
            continue
        pst   = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6, "label": particle})
        en    = summary["energy_MeV"].values
        res   = summary["ls_res_std"].values / summary["ls_res_mean"].values
        color = "#9467bd" if particle == "electron" else "#e377c2"
        ax.plot(en, res * 100, ls=pst["ls"], lw=pst["lw"], color=color,
                label=pst["label"])
    _add_energy_lines(ax)
    ax.set_xlabel("Initial energy (MeV)", fontsize=10)
    ax.set_ylabel("σ/μ  (%)", fontsize=10)
    ax.set_title("LS energy resolution\n(events reaching LS1)", fontsize=11)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    fig.suptitle("Liquid scintillator calorimetry", fontsize=13, y=1.01)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


def plot_containment(pdf: PdfPages, summaries: dict):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    ax = axes[0]
    for particle, summary in summaries.items():
        pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8, "label": particle})
        for col, color in [("trans_primInLS1", "#6a1d8a"),
                            ("trans_primInLS4", "#d7bde2")]:
            if col in summary.columns:
                layer = "LS1" if "LS1" in col else "LS4"
                ax.plot(summary["energy_MeV"].values, summary[col].values,
                        color=color, ls=pst["ls"], lw=pst["lw"],
                        label=f"{layer} ({pst['label']})")
    ax.set_ylim(-0.02, 1.05)
    ax.set_ylabel("Fraction of primaries", fontsize=10)
    _set_xaxis(ax)
    ax.set_title("LS layer entry fractions", fontsize=11)
    ax.legend(fontsize=8, ncol=2); ax.grid(True, alpha=0.3)

    ax = axes[1]
    for particle, summary in summaries.items():
        if "ls_total_mean" not in summary.columns:
            continue
        pst   = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8, "label": particle})
        en    = summary["energy_MeV"].values
        color = "#9467bd" if particle == "electron" else "#e377c2"
        ax.plot(en, summary["trans_primInLS1"].values if "trans_primInLS1" in summary.columns else np.nan,
                color=color, ls=pst["ls"], lw=pst["lw"], label=f"Enters LS1 ({pst['label']})")
        ax.plot(en, summary["ls_total_mean"].values / en,
                color=color, ls=":" , lw=1.2, label=f"LS edep/E₀ ({pst['label']})")
    ax.axhline(1.0, color="grey", lw=0.7, ls=":")
    ax.set_ylim(-0.05, 1.2)
    ax.set_ylabel("Fraction / ratio", fontsize=10)
    _set_xaxis(ax)
    ax.set_title("Containment quality", fontsize=11)
    ax.legend(fontsize=8, ncol=2); ax.grid(True, alpha=0.3)

    fig.suptitle("Containment analysis for LS calorimetry", fontsize=13, y=1.01)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


def plot_layer_profiles(pdf: PdfPages, all_groups: dict,
                        snapshot_energies: list):
    """Bar chart of edep per layer at snapshot energies, one panel per particle."""
    n_particles = len(all_groups)
    if n_particles == 0:
        return
    fig, axes = plt.subplots(1, n_particles,
                              figsize=(11 * n_particles / 2 + 3, 5),
                              sharey=True)
    if n_particles == 1:
        axes = [axes]

    layer_labels = [label for (_, _, label, _) in LAYERS]
    x = np.arange(len(LAYERS))
    colors_e = cm.viridis(np.linspace(0.2, 0.85, len(snapshot_energies)))

    for ax, (particle, groups) in zip(axes, all_groups.items()):
        snaps = find_snapshots(groups, snapshot_energies)
        bar_width = 0.7 / len(snaps)
        for i, (e_target, (e_actual, files)) in enumerate(sorted(snaps.items())):
            df = read_event_tree(files)
            if df.empty:
                continue
            vals = []
            for (edep_br, _, _, _) in LAYERS:
                vals.append(df[edep_br].mean() / 1e6 if edep_br in df.columns else 0.0)
            offset = (i - len(snaps) / 2 + 0.5) * bar_width
            ax.bar(x + offset, vals, bar_width * 0.9,
                   label=f"{e_actual:.3g} MeV", color=colors_e[i], alpha=0.85)
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(layer_labels, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("Mean edep (MeV/primary)", fontsize=10)
        ax.set_title(f"Layer profile — {PARTICLE_STYLES.get(particle,{}).get('label', particle)}",
                     fontsize=11)
        ax.legend(title="E₀", fontsize=8, title_fontsize=8)
        ax.grid(True, axis="y", which="both", alpha=0.3)

    fig.suptitle("Energy deposition profile at selected energies", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_angular_distributions(pdf: PdfPages, all_groups: dict,
                                snapshot_energies: list,
                                n_workers: int = None):
    """One row of histograms per particle. ClusterTree data loaded in parallel."""
    n_particles = len(all_groups)
    n_snaps     = len(snapshot_energies)
    colors      = cm.plasma(np.linspace(0.15, 0.85, n_snaps))

    # Build work list: (particle, e_actual, files) for all particles × snapshots
    work = []
    snap_map = {}   # particle → {e_target: (e_actual, files)}
    for particle, groups in all_groups.items():
        snaps = find_snapshots(groups, snapshot_energies)
        snap_map[particle] = snaps
        for e_target, (e_actual, files) in snaps.items():
            work.append((particle, e_actual, files))

    # Pre-fetch angles in parallel across all particles and snapshots
    print(f"  Loading angular data ({len(work)} cluster trees in parallel) ...")
    ang_cache = {}   # (particle, e_actual) → angles array
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=_MP_CTX) as ex:
        futures = {ex.submit(_angular_one, (e, f)): (p, e)
                   for p, e, f in work}
        for fut in as_completed(futures):
            p, e = futures[fut]
            try:
                _, angles = fut.result(timeout=300)   # 5-min cap per tree
            except Exception as err:
                print(f"  Warning angular {p} {e:.3g} MeV: {err}")
                angles = np.array([])
            ang_cache[(p, e)] = angles

    fig, axes = plt.subplots(n_particles, n_snaps,
                              figsize=(3.0 * n_snaps, 3.5 * n_particles),
                              sharey="row", sharex=False)
    if n_particles == 1:
        axes = axes[np.newaxis, :]
    if n_snaps == 1:
        axes = axes[:, np.newaxis]

    for row, (particle, groups) in enumerate(all_groups.items()):
        snaps = snap_map[particle]
        for col, ((e_target, (e_actual, _files)), color) in enumerate(
                zip(sorted(snaps.items()), colors)):
            ax     = axes[row, col]
            angles = ang_cache.get((particle, e_actual), np.array([]))
            if len(angles) == 0:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
                continue
            rms = np.sqrt(np.mean(angles**2))
            ax.hist(angles, bins=60, density=True, color=color, alpha=0.8)
            ax.axvline(rms, color="k", lw=1.1, ls="--", label=f"RMS={rms:.3f}°")
            ax.set_title(f"{e_actual:.3g} MeV", fontsize=9)
            ax.legend(fontsize=7)
            ax.set_xlim(left=0)
            if col == 0:
                plab = PARTICLE_STYLES.get(particle, {}).get("label", particle)
                ax.set_ylabel(f"{plab}\nDensity", fontsize=9)
            if row == n_particles - 1:
                ax.set_xlabel("Polar angle (degrees)", fontsize=8)

    fig.suptitle("Track angle in MM drift gas (multiple scattering spread)",
                 fontsize=11)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_angular_resolution(pdf: PdfPages, all_groups: dict, step: int = 10,
                             n_workers: int = None):
    """RMS polar angle vs energy — ClusterTree loaded in parallel across all particles."""
    fig, ax = plt.subplots(figsize=(8, 4.5))

    # Build work list across all particles at once for maximum parallelism
    work = []
    for particle, groups in all_groups.items():
        for energy in sorted(groups.keys())[::step]:
            work.append((particle, energy, groups[energy]))

    total_w = len(work)
    print(f"  Loading angular RMS data ({total_w} cluster trees in parallel) ...")
    ang_cache = {}   # (particle, energy) → angles
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=_MP_CTX) as ex:
        futures = {ex.submit(_angular_one, (e, f)): (p, e)
                   for p, e, f in work}
        done = 0
        for fut in as_completed(futures):
            done += 1
            p, e = futures[fut]
            try:
                _, angles = fut.result(timeout=300)   # 5-min cap per tree
            except Exception as err:
                print(f"  Warning angular RMS {p} {e:.4g} MeV: {err}")
                angles = np.array([])
            ang_cache[(p, e)] = angles
            if done % max(1, total_w // 5) == 0 or done == total_w:
                print(f"  ... angular {done}/{total_w} done")

    for particle, groups in all_groups.items():
        pst   = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6, "label": particle})
        color = "#1f77b4" if particle == "electron" else "#ff7f0e"
        en_vals, rms_vals = [], []
        for energy in sorted(groups.keys())[::step]:
            angles = ang_cache.get((particle, energy), np.array([]))
            if len(angles) < 20:
                continue
            rms_vals.append(np.sqrt(np.mean(angles**2)))
            en_vals.append(energy)
        if en_vals:
            ax.plot(en_vals, rms_vals, ls=pst["ls"], lw=pst["lw"],
                    color=color, marker="o", ms=3, label=pst["label"])

    ax.set_ylabel("RMS polar angle (degrees)", fontsize=11)
    _set_xaxis(ax)
    ax.set_title("Angular resolution from MS in upstream material\n"
                 "(straight-line fit in MM drift gas)", fontsize=11)
    ax.legend(fontsize=9); ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


# ── Detailed per-layer edep plot ─────────────────────────────────────────────

def plot_edep_detailed(pdf: PdfPages, summaries: dict):
    """
    Energy deposition in every individually scored material layer.
    Two panels stacked vertically so labels have room:
      Top : MM + PCB region (from Mylar to PCB Al foil)
      Bottom : scintillator + LS region
    Aggregation (e.g. total PCB) is done here in Python from the granular
    simulation branches.
    """
    mm_pcb_layers = [l for l in DETAILED_LAYERS
                     if l[0] not in ("edepHe3Gas", "edepLS1", "edepLS2",
                                     "edepLS3", "edepLS4", "edepLSCFRP",
                                     "edepScintWall", "edepScintTape", "edepScintAlFoil")]
    scint_ls_layers = [l for l in DETAILED_LAYERS
                       if l[0] in ("edepScintTape", "edepScintWall", "edepScintAlFoil",
                                   "edepLS1", "edepLS2", "edepLS3", "edepLS4", "edepLSCFRP")]

    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)

    for ax, layer_list, title in [
        (axes[0], mm_pcb_layers,  "Energy deposition: MM entrance + PCB"),
        (axes[1], scint_ls_layers,"Energy deposition: scintillator wall + LS stack"),
    ]:
        for (edep_br, _, label, color) in layer_list:
            col = f"edep_{edep_br}"
            for particle, summary in summaries.items():
                if col not in summary.columns:
                    continue
                vals = summary[col].values
                if np.nanmax(vals) < 1e-9:
                    continue
                pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6})
                ax.plot(summary["energy_MeV"].values, vals,
                        color=color, ls=pst["ls"], lw=pst["lw"])

        ax.set_yscale("log")
        ax.set_ylabel("Mean edep (MeV / primary)", fontsize=10)
        ax.set_title(title, fontsize=11)
        ax.grid(True, which="both", alpha=0.3)

        # Layer legend (color patches)
        active = [l for l in layer_list
                  if f"edep_{l[0]}" in next(iter(summaries.values())).columns
                  and next(iter(summaries.values()))[f"edep_{l[0]}"].max() > 1e-9]
        handles = [mpatches.Patch(color=l[3], label=l[2]) for l in active]
        leg1 = ax.legend(handles=handles, fontsize=7, ncol=2, loc="upper left",
                         title="Layer", title_fontsize=7)
        ax.add_artist(leg1)
        _particle_legend(ax, summaries, loc="lower right")

    _set_xaxis(axes[1])
    fig.suptitle("Granular energy deposition — all individually scored layers",
                 fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


# ── Coarse grouped plots ──────────────────────────────────────────────────────

def _coarse_legend(ax, summaries: dict, loc="lower right"):
    active = _active_trans_cols(summaries)
    handles = [mpatches.Patch(color=color, label=label)
               for (label, color, _, trans_col) in COARSE_GROUPS
               if f"trans_{trans_col}" in active]
    return ax.legend(handles=handles, loc=loc, fontsize=8,
                     title="Group", title_fontsize=8)


def plot_transmission_coarse(pdf: PdfPages, summaries: dict):
    active = _active_trans_cols(summaries)
    fig, ax = plt.subplots(figsize=(8, 5))
    for (label, color, _edep_brs, trans_col) in COARSE_GROUPS:
        col = f"trans_{trans_col}"
        if col not in active:
            continue
        for particle, summary in summaries.items():
            if col not in summary.columns:
                continue
            pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8})
            ax.plot(summary["energy_MeV"].values, summary[col].values,
                    color=color, ls=pst["ls"], lw=pst["lw"])
    ax.set_ylim(-0.02, 1.05)
    ax.axhline(1.0, color="grey", lw=0.6, ls=":")
    ax.set_ylabel("Transmission fraction", fontsize=11)
    _set_xaxis(ax)
    ax.set_title("Transmission by detector group", fontsize=12)
    ax.grid(True, alpha=0.3)
    leg1 = _coarse_legend(ax, summaries, loc="lower right")
    ax.add_artist(leg1)
    _particle_legend(ax, summaries, loc="center right")
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_edep_coarse(pdf: PdfPages, summaries: dict):
    active = _active_trans_cols(summaries)
    fig, ax = plt.subplots(figsize=(8, 5))
    for (label, color, edep_brs, trans_col) in COARSE_GROUPS:
        if f"trans_{trans_col}" not in active:
            continue
        for particle, summary in summaries.items():
            cols = [f"edep_{b}" for b in edep_brs if f"edep_{b}" in summary.columns]
            if not cols:
                continue
            total = sum(summary[c] for c in cols)
            if total.max() < 1e-8:
                continue
            pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8})
            ax.plot(summary["energy_MeV"].values, total.values,
                    color=color, ls=pst["ls"], lw=pst["lw"])
    ax.set_yscale("log")
    ax.set_ylabel("Mean energy deposition (MeV / primary)", fontsize=11)
    _set_xaxis(ax)
    ax.set_title("Energy deposition by detector group", fontsize=12)
    ax.grid(True, which="both", alpha=0.3)
    leg1 = _coarse_legend(ax, loc="upper left")
    ax.add_artist(leg1)
    _particle_legend(ax, summaries, loc="lower right")
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


# ── Survival vs path-length plots ─────────────────────────────────────────────

def _draw_stack_diagram(ax, xlim, skip_target: bool = False):
    """
    Thin schematic of the detector stack sharing the path-length x-axis.
    Air gaps take the figure background colour.  Very thin layers get a
    minimum display width.  Group labels with brackets appear above the blocks.

    skip_target: if True, the He-3 target blocks (z < 26.4 mm) are omitted,
                 e.g. for the sr90 calibration mode where the target is absent.
    """
    MIN_VIS_MM = 1.8

    ax.set_facecolor(ax.get_figure().get_facecolor())

    # ── Material blocks ───────────────────────────────────────────────────
    for (z0, z1, fc, tc, label) in DETECTOR_STACK:
        if skip_target and z1 <= 26.4:   # He-3 gas, Al wall, CFRP wall
            continue
        w_real = z1 - z0
        w_draw = max(w_real, MIN_VIS_MM)
        zc     = (z0 + z1) / 2
        rect   = mpatches.Rectangle(
            (zc - w_draw / 2, 0.04), w_draw, 0.72,
            facecolor=fc, edgecolor="white", linewidth=0.5, zorder=2)
        ax.add_patch(rect)

        if label and w_draw >= 3.0:
            rot = 0 if w_draw >= 14 else 90
            fs  = 7.0 if w_draw >= 14 else 5.5
            ax.text(zc, 0.40, label,
                    ha="center", va="center",
                    fontsize=fs, color=tc, rotation=rot,
                    fontweight="bold", zorder=3, clip_on=True)

    # ── Group labels with brackets above the blocks ───────────────────────
    # (z_start_mm, z_end_mm, label text)
    GROUPS = [
        (0.0,   26.4,  "Target"),
        (226.4, 262.6, "Micromegas"),
        (282.6, 285.8, "Scint.\nWall"),
        (305.8, 371.8, "Liquid Scint.\nCalorimeter"),
    ]
    bracket_y = 0.82   # top of bracket bar (in axes data coords)
    tick_h    = 0.07   # downward tick at each bracket end
    label_y   = 0.92   # text sits above the bracket bar

    for z0_g, z1_g, lbl in GROUPS:
        if skip_target and lbl == "Target":
            continue
        z0_c = max(z0_g, xlim[0])
        z1_c = min(z1_g, xlim[1])
        if z1_c <= z0_c:
            continue
        # horizontal bar
        ax.plot([z0_c, z1_c], [bracket_y, bracket_y],
                color="#444444", lw=1.1, clip_on=False, zorder=4)
        # end ticks pointing down toward the block
        for xv in (z0_c, z1_c):
            ax.plot([xv, xv], [bracket_y - tick_h, bracket_y],
                    color="#444444", lw=1.1, clip_on=False, zorder=4)
        # label
        ax.text((z0_c + z1_c) / 2, label_y, lbl,
                ha="center", va="bottom", fontsize=7.0,
                fontweight="bold", color="#222222",
                clip_on=False, zorder=5)

    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1.0)
    ax.set_autoscale_on(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(left=False, bottom=False,
                   labelleft=False, labelbottom=False)


def plot_survival(pdf: PdfPages, all_summaries: dict):
    """
    Two-page layout:
      Page 1 – detector diagram + full-width cumulative survival step function.
      Page 2 – detector diagram + full-width loss-rate (−dS/dz) plot.

    Both pages share the path-length x-axis with the diagram so that material
    blocks align with the data.  All particles are overlaid on each panel:
    solid = electron, dashed = positron.
    """
    XLIM        = (0.0, LS4_EXIT_MM + 6)
    DIAG_H      = 0.38   # height ratio: diagram row relative to data row
    snap_colors = cm.plasma(np.linspace(0.1, 0.85, len(SURVIVAL_ENERGIES_MEV)))

    # Detect whether the He-3 target was present in the simulation.
    # In sr90 mode trans_primInHe3Gas is always 0 — skip it to avoid a false
    # immediate drop to 0 in the survival curve.
    skip_target  = not _has_target(all_summaries)
    active_tcols = _active_trans_cols(all_summaries)

    # ── Pre-compute z/s arrays for every (particle, energy) combination ───
    all_curves = {}   # (particle, e_target) → (e_actual, z_pts, s_pts)
    for particle, summary in all_summaries.items():
        avail = summary["energy_MeV"].values
        for e_target in SURVIVAL_ENERGIES_MEV:
            idx      = np.argmin(np.abs(avail - e_target))
            e_actual = avail[idx]
            row_data = summary.iloc[idx]

            z_pts, s_pts = [0.0], [1.0]
            for (trans_col, z_entry, _) in SURVIVAL_COLS:
                col = f"trans_{trans_col}"
                if col not in active_tcols:
                    continue   # volume absent in this geometry
                if col in summary.columns:
                    val = row_data[col]
                    if not np.isnan(val):
                        z_pts.append(z_entry)
                        s_pts.append(float(val))

            all_curves[(particle, e_target)] = (
                e_actual, np.array(z_pts), np.array(s_pts))

    def _boundary_guides(ax):
        for (_, z_entry, _) in SURVIVAL_COLS:
            ax.axvline(z_entry, color="lightgrey", lw=0.7, ls=":", zorder=0)
        ax.axvline(LS4_EXIT_MM, color="lightgrey", lw=0.7, ls=":", zorder=0)

    def _make_legend(ax, particles):
        """Combined legend: colour = energy, linestyle = particle."""
        handles = []
        for color, e in zip(snap_colors, SURVIVAL_ENERGIES_MEV):
            handles.append(mlines.Line2D([], [], color=color, lw=2,
                                         label=f"{e:.3g} MeV"))
        handles.append(mlines.Line2D([], [], color="k", ls="-",  lw=1.5,
                                     label="e⁻  (solid)"))
        if "positron" in particles:
            handles.append(mlines.Line2D([], [], color="k", ls="--", lw=1.5,
                                         label="e⁺  (dashed)"))
        ax.legend(handles=handles, fontsize=7.5, ncol=2, loc="center left")

    # ── Page 1: cumulative survival ───────────────────────────────────────
    fig1 = plt.figure(figsize=(12, 5.5))
    gs1  = fig1.add_gridspec(2, 1, height_ratios=[DIAG_H, 1.0], hspace=0.05)

    ax_diag1 = fig1.add_subplot(gs1[0])
    _draw_stack_diagram(ax_diag1, xlim=XLIM, skip_target=skip_target)

    ax_s = fig1.add_subplot(gs1[1], sharex=ax_diag1)

    for particle, summary in all_summaries.items():
        pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8})
        for color, e_target in zip(snap_colors, SURVIVAL_ENERGIES_MEV):
            e_actual, z_pts, s_pts = all_curves[(particle, e_target)]
            # Extend to right edge so the last step is visible
            z_ext = np.append(z_pts, XLIM[1])
            s_ext = np.append(s_pts, s_pts[-1])
            ax_s.step(z_ext, s_ext, where="post",
                      color=color, ls=pst["ls"], lw=pst["lw"])

    _boundary_guides(ax_s)
    ax_s.set_ylim(-0.02, 1.06)
    ax_s.axhline(1.0, color="grey", lw=0.6, ls=":")
    ax_s.set_ylabel("Cumulative survival fraction", fontsize=11)
    ax_s.set_xlabel("Path from gun (mm)", fontsize=11)
    ax_s.grid(True, alpha=0.3)
    _make_legend(ax_s, list(all_summaries.keys()))

    fig1.suptitle("Primary particle survival through the detector stack",
                  fontsize=12, y=1.01)
    fig1.savefig(pdf, format="pdf", bbox_inches="tight")
    plt.close(fig1)

    # ── Page 2: loss rate ─────────────────────────────────────────────────
    fig2 = plt.figure(figsize=(12, 5.5))
    gs2  = fig2.add_gridspec(2, 1, height_ratios=[DIAG_H, 1.0], hspace=0.05)

    ax_diag2 = fig2.add_subplot(gs2[0])
    _draw_stack_diagram(ax_diag2, xlim=XLIM, skip_target=skip_target)

    ax_l = fig2.add_subplot(gs2[1], sharex=ax_diag2)

    for particle, summary in all_summaries.items():
        pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8})
        for color, e_target in zip(snap_colors, SURVIVAL_ENERGIES_MEV):
            e_actual, z_pts, s_pts = all_curves[(particle, e_target)]
            dz   = np.diff(z_pts)
            ds   = np.diff(s_pts)
            loss = -ds / dz
            z_mid = 0.5 * (z_pts[:-1] + z_pts[1:])
            ax_l.plot(z_mid, loss, color=color, ls=pst["ls"], lw=pst["lw"],
                      marker="o", ms=4)

    _boundary_guides(ax_l)
    ax_l.set_ylim(bottom=0)
    ax_l.set_ylabel("Loss rate  −dS/dz  (mm⁻¹)", fontsize=11)
    ax_l.set_xlabel("Path from gun (mm)", fontsize=11)
    ax_l.grid(True, alpha=0.3)
    _make_legend(ax_l, list(all_summaries.keys()))

    fig2.suptitle("Primary particle loss rate through the detector stack",
                  fontsize=12, y=1.01)
    fig2.savefig(pdf, format="pdf", bbox_inches="tight")
    plt.close(fig2)


# ── Trigger and calorimetry plots ────────────────────────────────────────────

def _plateau_mean(summary: pd.DataFrame, col: str, top_frac: float = 0.15) -> float:
    """Mean of the highest top_frac of values — used as plateau estimate."""
    vals = summary[col].dropna().values
    n_hi = max(1, int(len(vals) * top_frac))
    return float(np.mean(np.sort(vals)[-n_hi:]))


def compute_trigger_fractions(groups: dict,
                              thresh_scint_eV: float,
                              thresh_ls1_eV: float,
                              n_workers: int = None) -> dict:
    """
    Per-energy coincidence trigger fraction — parallelised over energy points.
    Reads only two EventTree branches for speed.
    """
    items = [(e, f, thresh_scint_eV, thresh_ls1_eV) for e, f in groups.items()]
    n = len(items)
    result = {}
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=_MP_CTX) as ex:
        futures = {ex.submit(_trigger_one, a): a[0] for a in items}
        done = 0
        for fut in as_completed(futures):
            done += 1
            energy = futures[fut]
            try:
                e, frac = fut.result()
                result[e] = frac
            except Exception as err:
                print(f"  Warning trigger {energy:.4g} MeV: {err}")
                result[energy] = np.nan
            if done % max(1, n // 10) == 0 or done == n:
                print(f"  ... trigger {done}/{n} done")
    return result


def plot_trigger(pdf: PdfPages, all_summaries: dict, all_groups: dict,
                 plateau_frac: float = 0.30, n_workers: int = None):
    """
    Left panel : mean edep in plastic scint. and LS1 vs energy.
                 Horizontal dashed lines mark the coincidence thresholds
                 (plateau_frac × high-energy plateau for each detector).
    Right panel: coincidence trigger fraction vs energy.
    Both electrons and positrons overlaid.
    """
    fig, (ax_e, ax_t) = plt.subplots(1, 2, figsize=(12, 4.5))

    thresholds = {}   # particle → (thresh_scint_eV, thresh_ls1_eV)

    for particle, summary in all_summaries.items():
        pst   = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8, "label": particle})
        en    = summary["energy_MeV"].values

        col_s = "edep_edepScintWall"
        col_l = "edep_edepLS1"
        if col_s not in summary.columns or col_l not in summary.columns:
            continue

        vals_s = summary[col_s].values   # MeV
        vals_l = summary[col_l].values   # MeV

        ax_e.plot(en, vals_s, color="#d62728", ls=pst["ls"], lw=pst["lw"])
        ax_e.plot(en, vals_l, color="#6a1d8a", ls=pst["ls"], lw=pst["lw"])

        # Thresholds from plateau
        plat_s = _plateau_mean(summary, col_s)   # MeV
        plat_l = _plateau_mean(summary, col_l)   # MeV
        thr_s  = plateau_frac * plat_s           # MeV
        thr_l  = plateau_frac * plat_l           # MeV
        thresholds[particle] = (thr_s * 1e6, thr_l * 1e6)   # eV for uproot

    # Draw threshold lines once (same for both particles)
    if thresholds:
        thr_s_MeV = plateau_frac * _plateau_mean(
            next(iter(all_summaries.values())), "edep_edepScintWall")
        thr_l_MeV = plateau_frac * _plateau_mean(
            next(iter(all_summaries.values())), "edep_edepLS1")
        ax_e.axhline(thr_s_MeV, color="#d62728", lw=1.0, ls=":",
                     label=f"Thr. scint. {thr_s_MeV*1000:.0f} keV")
        ax_e.axhline(thr_l_MeV, color="#6a1d8a", lw=1.0, ls=":",
                     label=f"Thr. LS1 {thr_l_MeV*1000:.0f} keV")

    ax_e.set_yscale("log")
    ax_e.set_ylabel("Mean energy deposition (MeV)", fontsize=10)
    _set_xaxis(ax_e)
    ax_e.set_title("Detector response: plastic scint. (red) and LS1 (purple)",
                   fontsize=10)
    ax_e.legend(fontsize=8)
    ax_e.grid(True, which="both", alpha=0.3)

    # Layer color patches for legend
    handles_det = [mpatches.Patch(color="#d62728", label="Plastic scint."),
                   mpatches.Patch(color="#6a1d8a", label="Liq. scint. 1")]
    leg1 = ax_e.legend(handles=handles_det, loc="upper left", fontsize=8)
    ax_e.add_artist(leg1)
    # Threshold lines legend
    ax_e.legend([ax_e.lines[-2], ax_e.lines[-1]],
                [ax_e.lines[-2].get_label(), ax_e.lines[-1].get_label()],
                fontsize=7, loc="lower right")

    # Trigger fractions
    for particle, (thr_s_eV, thr_l_eV) in thresholds.items():
        if particle not in all_groups:
            continue
        pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8, "label": particle})
        print(f"  Computing trigger fractions for {particle} "
              f"(scint>{thr_s_eV/1e6*1000:.0f} keV, LS1>{thr_l_eV/1e6*1000:.0f} keV)...")
        frac = compute_trigger_fractions(all_groups[particle], thr_s_eV, thr_l_eV,
                                         n_workers=n_workers)
        en_t = np.array(sorted(frac.keys()))
        fr_t = np.array([frac[e] for e in en_t])
        ax_t.plot(en_t, fr_t, color="#9467bd", ls=pst["ls"], lw=pst["lw"],
                  label=pst["label"])

    ax_t.set_ylim(-0.02, 1.05)
    ax_t.axhline(1.0, color="grey", lw=0.6, ls=":")
    ax_t.set_ylabel("Coincidence trigger fraction", fontsize=10)
    _set_xaxis(ax_t)
    ax_t.set_title(f"Trigger: plastic scint. AND LS1 above threshold\n"
                   f"(threshold = {plateau_frac*100:.0f}% of high-energy plateau)",
                   fontsize=10)
    _particle_legend(ax_t, all_summaries, loc="lower right")
    ax_t.grid(True, alpha=0.3)

    fig.suptitle("Coincidence trigger: plastic scintillator ∩ LS1", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_calorimetry_fractions(pdf: PdfPages, summaries: dict):
    """
    Fraction of primaries in three mutually exclusive categories per energy:
      • No entry to LS   = 1 − primInLS1           (missed calorimeter)
      • Contained in LS  = primInLS1 − primInLSCFRP5 (good calorimetry)
      • Punch-through    = primInLSCFRP5             (exited back wall)
    If primInLSCFRP5 is unavailable the punch-through series is omitted.
    """
    fig, ax = plt.subplots(figsize=(9, 5))

    for particle, summary in summaries.items():
        pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8, "label": particle})
        en  = summary["energy_MeV"].values

        if "trans_primInLS1" not in summary.columns:
            continue
        frac_in   = summary["trans_primInLS1"].values
        frac_miss = 1.0 - frac_in
        have_pt   = "trans_primInLSCFRP5" in summary.columns
        frac_pt   = summary["trans_primInLSCFRP5"].values if have_pt else np.zeros_like(frac_in)
        frac_cont = frac_in - frac_pt

        lbl = pst["label"]
        ax.plot(en, frac_miss, color="#d62728", ls=pst["ls"], lw=pst["lw"],
                label=f"No entry to LS ({lbl})")
        ax.plot(en, frac_cont, color="#2ca02c", ls=pst["ls"], lw=pst["lw"],
                label=f"Contained ({lbl})")
        if have_pt:
            ax.plot(en, frac_pt, color="#ff7f0e", ls=pst["ls"], lw=pst["lw"],
                    label=f"Punch-through ({lbl})")

    ax.set_ylim(-0.02, 1.05)
    ax.axhline(1.0, color="grey", lw=0.6, ls=":")
    ax.set_ylabel("Fraction of primaries", fontsize=11)
    _set_xaxis(ax)
    ax.set_title("Calorimetry containment fractions\n"
                 "Contained = entered LS1 and stopped before back CFRP wall",
                 fontsize=11)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_linearity_standalone(pdf: PdfPages, summaries: dict):
    """
    Standalone energy linearity plot: total LS active edep vs initial energy.
    Linear axes. Both particles. Includes CFRP dead-layer losses as annotation.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax_lin = axes[0]   # main linearity
    ax_res = axes[1]   # residual / fraction

    e_max = 0.0
    for particle, summary in summaries.items():
        if "ls_total_mean" not in summary.columns:
            continue
        pst   = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8, "label": particle})
        color = "#9467bd" if particle == "electron" else "#e377c2"
        en    = summary["energy_MeV"].values
        mean  = summary["ls_total_mean"].values
        std   = summary["ls_total_std"].values
        e_max = max(e_max, en.max())

        # Shaded ±1σ band
        ax_lin.fill_between(en, mean - std, mean + std,
                             alpha=0.15, color=color)
        ax_lin.plot(en, mean, color=color, ls=pst["ls"], lw=pst["lw"],
                    label=pst["label"])

        # Residual: (edep_LS - E₀) / E₀
        with np.errstate(divide="ignore", invalid="ignore"):
            res = (mean - en) / en
        ax_res.plot(en, res, color=color, ls=pst["ls"], lw=pst["lw"],
                    label=pst["label"])

        # CFRP dead losses if available
        if "edep_edepLSCFRP" in summary.columns:
            cfrp = summary["edep_edepLSCFRP"].values
            ax_lin.plot(en, mean + cfrp, color=color, ls=":",
                        lw=1.0, alpha=0.6, label=f"{pst['label']} (+ CFRP losses)")

    # Ideal response
    e_ref = np.linspace(0, e_max, 200)
    ax_lin.plot(e_ref, e_ref, "k--", lw=0.9, label="Ideal (E_dep = E₀)")
    ax_lin.set_xlim(0, e_max * 1.02)
    ax_lin.set_ylim(0, e_max * 1.02)
    ax_lin.set_xlabel("Initial energy (MeV)", fontsize=11)
    ax_lin.set_ylabel("Total LS active edep (MeV)", fontsize=11)
    ax_lin.set_title("Calorimeter response — linear scale", fontsize=12)
    _add_energy_lines(ax_lin)
    ax_lin.legend(fontsize=8)
    ax_lin.grid(True, alpha=0.3)

    ax_res.axhline(0, color="k", lw=0.8, ls="--")
    ax_res.set_xlabel("Initial energy (MeV)", fontsize=11)
    ax_res.set_ylabel("(E_dep − E₀) / E₀", fontsize=11)
    ax_res.set_title("Relative deviation from ideal response", fontsize=12)
    _add_energy_lines(ax_res)
    ax_res.legend(fontsize=8)
    ax_res.grid(True, alpha=0.3)

    fig.suptitle("LS calorimeter energy linearity\n"
                 "Shaded band = ±1σ event-to-event spread", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Full-experiment stack analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--indir",    required=True)
    p.add_argument("--gas",      default="ArIso")
    p.add_argument("--particles", nargs="+", default=["electron", "positron"],
                   help="Particle types to load (space-separated)")
    p.add_argument("--outfile",  default="full_experiment_analysis.pdf")
    p.add_argument("--angular-step", type=int, default=10,
                   help="Process every Nth energy for angular RMS plot")
    p.add_argument("--no-angular", action="store_true",
                   help="Skip angular analysis (faster)")
    p.add_argument("--workers", type=int, default=None,
                   help="Parallel worker processes (default: all CPU cores)")
    p.add_argument("--prefix", default="full",
                   help="File prefix: 'full' for full-experiment run, "
                        "'sr90' for Sr-90 calibration run")
    return p.parse_args()


def main():
    args = parse_args()

    plt.rcParams.update({
        "figure.dpi": 150, "axes.titlesize": 11,
        "axes.labelsize": 10, "xtick.labelsize": 9,
        "ytick.labelsize": 9, "legend.fontsize": 9,
    })

    nw = args.workers   # None → ProcessPoolExecutor uses os.cpu_count()
    print(f"Workers: {nw or os.cpu_count()} (set with --workers)")

    all_groups    = {}
    all_summaries = {}

    for particle in args.particles:
        print(f"\nLoading {particle} ...")
        groups = group_files(args.indir, args.gas, particle, prefix=args.prefix)
        if not groups:
            continue
        print(f"  {len(groups)} energy points: "
              f"{min(groups):.4g} – {max(groups):.4g} MeV")
        all_groups[particle] = groups
        print(f"  Computing summaries ...")
        all_summaries[particle] = compute_summary(groups, n_workers=nw)

    if not all_summaries:
        sys.exit("No data loaded.")

    # Save summary CSVs
    for particle, summary in all_summaries.items():
        csv = Path(args.outfile).with_name(
            Path(args.outfile).stem + f"_{particle}.csv")
        summary.to_csv(csv, index=False)
        print(f"Summary saved → {csv}")

    e_range = {p: f"{min(g):.4g}–{max(g):.4g} MeV"
               for p, g in all_groups.items()}

    print(f"\nWriting {args.outfile} ...")
    with PdfPages(args.outfile) as pdf:

        # Cover
        fig = plt.figure(figsize=(8, 4))
        fig.text(0.5, 0.68, "Full-experiment stack analysis",
                 ha="center", fontsize=18, fontweight="bold")
        fig.text(0.5, 0.52,
                 f"Gas: {args.gas}\n" +
                 "\n".join(f"{p}: {r}" for p, r in e_range.items()) +
                 f"\nInput: {args.indir}",
                 ha="center", fontsize=10, color="grey")
        pdf.savefig(fig); plt.close(fig)

        print("  Transmission (all layers) ...")
        plot_transmission(pdf, all_summaries)
        print("  Transmission (simplified: drift / plastic scint. / LS1 / LS4) ...")
        plot_transmission_simple(pdf, all_summaries)
        print("  Transmission (coarse: detector groups) ...")
        plot_transmission_coarse(pdf, all_summaries)

        print("  Edep (all layers) ...")
        plot_edep(pdf, all_summaries)
        print("  Edep (simplified) ...")
        plot_edep_simple(pdf, all_summaries)
        print("  Edep (coarse: detector groups) ...")
        plot_edep_coarse(pdf, all_summaries)
        print("  Edep (granular: every individual layer) ...")
        plot_edep_detailed(pdf, all_summaries)

        print("  Trigger: plastic scint. + LS1 coincidence ...")
        plot_trigger(pdf, all_summaries, all_groups, n_workers=nw)

        print("  Calorimetry fractions (no-entry / contained / punch-through) ...")
        plot_calorimetry_fractions(pdf, all_summaries)

        print("  Calorimetry linearity (standalone, linear axes) ...")
        plot_linearity_standalone(pdf, all_summaries)

        print("  Calorimetry (full 3-panel) ...")
        plot_calorimetry(pdf, all_summaries)

        print("  Containment ...")
        plot_containment(pdf, all_summaries)

        print("  Survival vs path length ...")
        plot_survival(pdf, all_summaries)

        print("  Layer profiles ...")
        plot_layer_profiles(pdf, all_groups, SNAPSHOT_ENERGIES)

        if not args.no_angular:
            print("  Angular distributions ...")
            plot_angular_distributions(pdf, all_groups, SNAPSHOT_ENERGIES,
                                       n_workers=nw)
            print(f"  Angular RMS (every {args.angular_step}th point) ...")
            plot_angular_resolution(pdf, all_groups, step=args.angular_step,
                                    n_workers=nw)
        else:
            print("  Skipping angular analysis (--no-angular)")

    print(f"\nDone → {args.outfile}")


if __name__ == "__main__":
    # The __main__ guard is required for ProcessPoolExecutor on macOS/Windows.
    # On Linux (lxplus) it is not strictly necessary but is good practice.
    main()
