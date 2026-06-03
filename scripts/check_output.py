#!/usr/bin/env python3
"""
check_output.py — Quick sanity check for a single mx17_full_sim ROOT output file.

Usage:
    python check_output.py path/to/x17_ipc_pairs_jobXXX_t0.root [-o report.pdf]

Prints a text report and (optionally) saves a 1-page PDF of key distributions.
Pass/fail checks are printed at the end so issues are immediately obvious.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import uproot
import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ── helpers ────────────────────────────────────────────────────────────────────

def _decode(series):
    if series.empty:
        return series
    s = series.iloc[0]
    if isinstance(s, bytes):
        return series.apply(lambda x: x.rstrip(b"\x00").decode("ascii", errors="ignore"))
    if hasattr(series, "str"):
        return series.str.rstrip("\x00")
    return series.astype(str)


def banner(title):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Sanity check a single mx17_full_sim ROOT file")
    ap.add_argument("rootfile", help="Path to ROOT file (e.g. x17_ipc_pairs_job096_t0.root)")
    ap.add_argument("-o", "--output", default=None,
                    help="Save a PDF of plots (optional)")
    args = ap.parse_args()

    path = Path(args.rootfile)
    if not path.exists():
        print(f"ERROR: file not found: {path}")
        sys.exit(1)

    print(f"\nFile: {path}  ({path.stat().st_size / 1e6:.1f} MB)")

    f = uproot.open(str(path))

    passes = []   # (label, bool, details)

    # ── 1. Tree inventory ──────────────────────────────────────────────────────
    banner("Tree inventory")
    all_keys = f.keys(cycle=True)
    tree_info = {}
    for k in all_keys:
        obj = f[k]
        n = getattr(obj, "num_entries", None)
        name, cycle = k.rsplit(";", 1) if ";" in k else (k, "?")
        print(f"  {k:<30}  entries={n}")
        if name not in tree_info:
            tree_info[name] = {"cycles": [], "entries": []}
        tree_info[name]["cycles"].append(int(cycle) if cycle.isdigit() else cycle)
        tree_info[name]["entries"].append(n)

    has_hit  = "HitTree"   in tree_info
    has_evt  = "EventTree" in tree_info
    passes.append(("HitTree present",   has_hit,  ""))
    passes.append(("EventTree present", has_evt,  ""))

    # Warn about multiple cycles (indicates accumulated runs)
    for name, info in tree_info.items():
        if len(info["cycles"]) > 1:
            print(f"\n  WARNING: {name} has {len(info['cycles'])} cycles "
                  f"(runs accumulated without RECREATE). "
                  f"Latest cycle will be used.")
            passes.append((f"{name} single cycle", False,
                           f"{len(info['cycles'])} cycles found — file was written multiple times"))

    # ── 2. EventTree check ─────────────────────────────────────────────────────
    banner("EventTree")
    evt_df = None
    if has_evt:
        et = f["EventTree"]
        print(f"  Entries : {et.num_entries}")
        print(f"  Branches: {et.keys()}")
        evt_df = et.arrays(library="pd")

        print(f"\n  First 5 rows:")
        print(evt_df.head(5).to_string(index=False))

        # Stats on numeric columns
        print(f"\n  Value ranges:")
        for col in evt_df.columns:
            mn, mx = evt_df[col].min(), evt_df[col].max()
            print(f"    {col:<20} min={mn:>12.4g}   max={mx:>12.4g}")

        # Pass/fail: non-zero physics
        all_zero = (evt_df.drop(columns=["eventID"]).abs().max() == 0).all()
        passes.append(("EventTree physics non-zero", not all_zero,
                       "ALL physics branches are exactly 0 — EventData was never filled. "
                       "Check binary version or RECREATE mode."))

        # event_type distribution
        if "event_type" in evt_df.columns:
            vc = evt_df["event_type"].value_counts().sort_index().to_dict()
            print(f"\n  event_type distribution: {vc}")
            has_ipc = 1 in vc and vc[1] > 0
            passes.append(("IPC events present (event_type=1)", has_ipc,
                           f"All events are X17 (event_type=0 only): {vc}"))
            frac_ipc = vc.get(1, 0) / et.num_entries if et.num_entries else 0
            print(f"  IPC fraction: {frac_ipc:.3f} (expected ~0.5 for default --ipc 0.5)")

        # eventID continuity
        if "eventID" in evt_df.columns:
            eids = evt_df["eventID"].values
            expected = np.arange(len(eids))
            ids_ok = np.array_equal(np.sort(eids), expected)
            passes.append(("EventTree eventIDs sequential", ids_ok,
                           f"eventID range: {eids.min()} to {eids.max()}"))

        # inv_mass non-zero
        if "inv_mass" in evt_df.columns:
            n_zero_mass = (evt_df["inv_mass"] == 0).sum()
            passes.append(("inv_mass non-zero", n_zero_mass == 0,
                           f"{n_zero_mass}/{len(evt_df)} events have inv_mass=0"))
    else:
        print("  (not found)")

    # ── 3. HitTree check ──────────────────────────────────────────────────────
    banner("HitTree")
    hit_df = None
    if has_hit:
        ht = f["HitTree"]
        print(f"  Entries : {ht.num_entries}")
        print(f"  Branches: {ht.keys()}")

        # Load a sample
        nread = min(ht.num_entries, 500_000)
        hit_df = ht.arrays(["eventID","detType","particle","ke","edep","armID"],
                           entry_stop=nread, library="pd")

        # Decode string columns
        hit_df["detType"]  = _decode(hit_df["detType"])
        hit_df["particle"] = _decode(hit_df["particle"])

        # detType distribution
        det_counts = hit_df["detType"].value_counts()
        print(f"\n  detType counts (first {nread//1000}k hits):")
        for det, cnt in det_counts.items():
            print(f"    {det:<20} {cnt:>8,}")

        expected_dets = {"DriftGas", "PlasticScint", "LiqScint_1", "LiqScint_2"}
        found_dets    = set(det_counts.index)
        missing_dets  = expected_dets - found_dets
        passes.append(("Key detTypes present", len(missing_dets) == 0,
                       f"Missing: {missing_dets}" if missing_dets else ""))

        # particle distribution
        part_counts = hit_df["particle"].value_counts()
        print(f"\n  Particle counts:")
        for p, cnt in part_counts.items():
            print(f"    {p:<20} {cnt:>8,}")

        has_ep = "e+" in part_counts.index
        has_em = "e-" in part_counts.index
        passes.append(("e+ and e- hits present", has_ep and has_em,
                       f"particles found: {list(part_counts.index[:8])}"))

        # KE range for e+/e-
        primary = hit_df[hit_df["particle"].isin(["e+","e-"])]
        if not primary.empty:
            print(f"\n  e+/e- KE range:")
            for p in ["e+","e-"]:
                sub = primary[primary["particle"]==p]["ke"]
                if not sub.empty:
                    print(f"    {p}: [{sub.min():.3f}, {sub.max():.3f}] MeV  "
                          f"mean={sub.mean():.3f} MeV")

        # arm distribution
        if "armID" in hit_df.columns:
            arm_counts = hit_df["armID"].value_counts().sort_index()
            print(f"\n  Hits per arm:")
            for arm, cnt in arm_counts.items():
                print(f"    arm {arm}: {cnt:,}")
            n_arms = len(arm_counts)
            passes.append(("All 4 arms have hits", n_arms == 4,
                           f"Only {n_arms} arms have hits: {list(arm_counts.index)}"))

        # edep sanity
        neg_edep = (hit_df["edep"] < 0).sum()
        passes.append(("No negative edep", neg_edep == 0,
                       f"{neg_edep} hits with edep<0"))

        # HitTree eventID range
        eid_min = hit_df["eventID"].min()
        eid_max = hit_df["eventID"].max()
        print(f"\n  eventID range in sample: {eid_min} to {eid_max}")
        # Cross-check with EventTree
        if evt_df is not None:
            evt_ids = set(evt_df["eventID"].values)
            hit_ids = set(hit_df["eventID"].unique())
            orphan_hit_events = hit_ids - evt_ids
            if orphan_hit_events:
                passes.append(("HitTree eventIDs subset of EventTree", False,
                               f"{len(orphan_hit_events)} HitTree eventIDs not in EventTree "
                               f"(sample: {sorted(orphan_hit_events)[:5]})"))
            else:
                passes.append(("HitTree eventIDs subset of EventTree", True, ""))
    else:
        print("  (not found)")

    # ── 4. Pass/fail summary ───────────────────────────────────────────────────
    banner("Pass / Fail Summary")
    ok_count  = sum(1 for _, ok, _ in passes if ok)
    fail_count = sum(1 for _, ok, _ in passes if not ok)
    for label, ok, detail in passes:
        sym = "PASS" if ok else "FAIL"
        line = f"  [{sym}]  {label}"
        if not ok and detail:
            line += f"\n         ↳ {detail}"
        print(line)
    print(f"\n  {ok_count} passed, {fail_count} failed")

    # ── 5. Plots (optional) ────────────────────────────────────────────────────
    if args.output and HAS_MPL and hit_df is not None:
        out_path = Path(args.output)
        print(f"\nSaving plots to {out_path} ...")
        with PdfPages(str(out_path)) as pdf:
            fig, axes = plt.subplots(2, 3, figsize=(15, 9))
            fig.suptitle(f"Sanity check: {path.name}", fontsize=11)

            ax = axes[0, 0]
            if evt_df is not None and "event_type" in evt_df.columns and not all_zero:
                vc = evt_df["event_type"].value_counts().sort_index()
                ax.bar([str(k) for k in vc.index], vc.values,
                       color=["#e84040","#4a90d9","#888888"])
                ax.set_title("event_type distribution")
                ax.set_xlabel("event_type (0=X17, 1=IPC, -1=single)")
            else:
                ax.text(0.5, 0.5, "EventTree empty\n(all zeros)",
                        ha="center", va="center", color="red", fontsize=14,
                        transform=ax.transAxes)
                ax.set_title("event_type distribution")

            ax = axes[0, 1]
            if evt_df is not None and "em_ke" in evt_df.columns and not all_zero:
                for col, label, color in [("em_ke","e- truth KE","#4a90d9"),
                                           ("ep_ke","e+ truth KE","#e84040")]:
                    vals = evt_df[col].values
                    vals = vals[vals > 0]
                    if len(vals):
                        ax.hist(vals, bins=80, histtype="step", label=label, color=color)
                ax.set_xlabel("Kinetic energy [MeV]")
                ax.set_title("Truth KE from EventTree")
                ax.legend(fontsize=8)
            else:
                ax.text(0.5, 0.5, "No EventTree data", ha="center", va="center",
                        transform=ax.transAxes)
                ax.set_title("Truth KE (not available)")

            ax = axes[0, 2]
            if evt_df is not None and "openingAngle" in evt_df.columns and not all_zero:
                angles = evt_df["openingAngle"].values
                angles = angles[angles > 0]
                ax.hist(angles, bins=80, color="#2ca02c", histtype="step")
                ax.set_xlabel("Opening angle [deg]")
                ax.set_title("Truth opening angle")
            else:
                ax.text(0.5, 0.5, "No EventTree data", ha="center", va="center",
                        transform=ax.transAxes)
                ax.set_title("Opening angle (not available)")

            ax = axes[1, 0]
            det_order = ["DriftGas","AmpGas","PlasticScint",
                         "LiqScint_1","LiqScint_2","BackScintL","BackScintR"]
            det_counts_all = hit_df["detType"].value_counts()
            det_counts_plot = {d: det_counts_all.get(d, 0) for d in det_order
                               if d in det_counts_all.index}
            if det_counts_plot:
                ax.barh(list(det_counts_plot.keys()),
                        [v/1e3 for v in det_counts_plot.values()],
                        color="#ff7f0e")
                ax.set_xlabel("Hits (×10³)")
                ax.set_title(f"Hits per detType (first {nread//1000}k rows)")

            ax = axes[1, 1]
            primary2 = hit_df[hit_df["particle"].isin(["e+","e-"])]
            if not primary2.empty:
                for p, color in [("e+","#e84040"), ("e-","#4a90d9")]:
                    sub = primary2[primary2["particle"]==p]["ke"]
                    if not sub.empty:
                        ax.hist(sub.values, bins=100, histtype="step",
                                label=p, color=color, density=True)
                ax.set_xlabel("KE at hit [MeV]")
                ax.set_title("e±  KE spectrum in HitTree")
                ax.legend(fontsize=8)

            ax = axes[1, 2]
            if "armID" in hit_df.columns:
                arm_counts2 = hit_df["armID"].value_counts().sort_index()
                ax.bar([f"arm {i}" for i in arm_counts2.index],
                       arm_counts2.values / 1e3, color="#9467bd")
                ax.set_ylabel("Hits (×10³)")
                ax.set_title("Hits per arm")

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            # Second page: edep distributions per layer
            fig, axes = plt.subplots(2, 3, figsize=(15, 8))
            fig.suptitle("Energy deposition by layer", fontsize=11)
            det_colors = {
                "DriftGas":    "#2ca02c",
                "PlasticScint":"#ff7f0e",
                "LiqScint_1":  "#6a1d8a",
                "LiqScint_2":  "#9b59b6",
                "BackScintL":  "#d62728",
                "BackScintR":  "#e57373",
            }
            plot_dets = [d for d in det_order if d in det_counts_all.index and d != "AmpGas"]
            for i, det in enumerate(plot_dets[:6]):
                ax = axes[i // 3][i % 3]
                sub = hit_df[hit_df["detType"] == det]["edep"]
                sub = sub[sub > 0]
                if not sub.empty:
                    ax.hist(sub.values, bins=80, histtype="stepfilled",
                            color=det_colors.get(det, "gray"), alpha=0.7)
                    ax.set_xlabel("edep [eV]")
                    ax.set_yscale("log")
                ax.set_title(det)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        print(f"Saved: {out_path}")
    elif args.output and not HAS_MPL:
        print("WARNING: matplotlib not available, skipping PDF output")

    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
