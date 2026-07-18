#!/usr/bin/env python3
"""
make_capture_library.py — capture-vertex library + budget from a neutron run
=============================================================================
Reads EventTree from `mx17_full_sim --neutron` output files and writes:

  1. A CSV capture-vertex library for the biased gamma-source mode
     (`mx17_full_sim --gamma-source <csv>`): one row per Al/CFRP capture,
     columns volume,x_mm,y_mm,z_mm.

  2. A capture-budget summary (stdout): captures per volume per incident
     neutron, split by decade of neutron energy — the normalisation input
     for the gamma-source run weights.

Usage:
    python3 make_capture_library.py /eos/.../neutrons/*.root \
        -o capture_lib.csv [--max-rows 500000]
"""

import argparse
import glob
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import uproot

WALL_VOLS = ("He3Cap_Al", "He3Cap_CFRP")


def collect_files(inputs):
    files = []
    for inp in inputs:
        p = Path(inp)
        if p.is_dir():
            files.extend(sorted(p.glob("*.root")))
        elif "*" in inp:
            files.extend(sorted(glob.glob(inp)))
        else:
            files.append(p)
    return [str(f) for f in files if Path(f).exists()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("-o", "--outfile", default="capture_lib.csv")
    ap.add_argument("--max-rows", type=int, default=500_000,
                    help="cap on library rows (uniform thinning above)")
    ap.add_argument("--gas-lib", default=None, metavar="GAS_CSV",
                    help="also write a He3Gas capture-position library "
                         "(pair-vertex sampling: mx17_full_sim --pair-vertex-lib)")
    ap.add_argument("--gas-emax", type=float, default=2.0,
                    help="max neutron energy [eV] for --gas-lib rows "
                         "(default 2 eV ~ TOF > 1 ms at 19.5 m)")
    args = ap.parse_args()

    files = collect_files(args.inputs)
    if not files:
        print("No input files found", file=sys.stderr)
        sys.exit(1)

    n_total = 0
    vol_counts = Counter()
    # capture counts per (volume, log10(E) decade)
    decade_counts = Counter()
    lib_rows = []   # (vol, x, y, z)
    gas_rows = []   # He3Gas capture positions with E_n <= gas_emax

    def _clean(arr):
        return np.array([v.rstrip("\x00") if isinstance(v, str)
                         else v.rstrip(b"\x00").decode() for v in arr])

    for fp in files:
        with uproot.open(fp) as f:
            if "EventTree" not in f:
                continue
            t = f["EventTree"].arrays(
                ["event_type", "neutron_E_eV", "capture_vol", "capture_proc",
                 "cap_x", "cap_y", "cap_z"], library="np")
        is_n = t["event_type"] == 2
        n_total += int(is_n.sum())
        vols  = _clean(t["capture_vol"][is_n])
        procs = _clean(t["capture_proc"][is_n])
        e_ev = t["neutron_E_eV"][is_n]
        with np.errstate(divide="ignore"):
            dec = np.floor(np.log10(np.clip(e_ev, 1e-30, None))).astype(int)

        for v, p, d in zip(vols, procs, dec):
            if v:
                vol_counts[f"{v} [{p}]"] += 1
                decade_counts[(v, int(d))] += 1

        # Library: RADIATIVE captures in the walls (those emit the γ cascade)
        wall = np.isin(vols, WALL_VOLS) & (procs == "nCapture")
        for v, x, y, z in zip(vols[wall], t["cap_x"][is_n][wall],
                              t["cap_y"][is_n][wall], t["cap_z"][is_n][wall]):
            lib_rows.append((v, x, y, z))

        # Gas library: any terminal capture in He3Gas below --gas-emax.  The
        # (n,p) positions trace the self-shielded absorption profile, which is
        # the same spatial law the rare (n,γ) follows — so they serve as the
        # pair-vertex distribution for thermal X17/IPC generation.
        if args.gas_lib:
            gas = (vols == "He3Gas") & (e_ev <= args.gas_emax)
            for x, y, z in zip(t["cap_x"][is_n][gas], t["cap_y"][is_n][gas],
                               t["cap_z"][is_n][gas]):
                gas_rows.append(("He3Gas", x, y, z))

    if n_total == 0:
        print("No neutron events (event_type==2) found", file=sys.stderr)
        sys.exit(1)

    print(f"Neutrons processed : {n_total:,}")
    print(f"\nCapture budget (per incident neutron):")
    for v, c in vol_counts.most_common():
        print(f"  {v:<20s}: {c:>10,}   ({c/n_total:.3e}/n)")
    n_cap = sum(vol_counts.values())
    print(f"  {'no capture/escape':<20s}: {n_total-n_cap:>10,}")

    print(f"\nWall captures by neutron-energy decade (log10 E/eV):")
    for v in WALL_VOLS:
        decs = sorted(d for (vv, d) in decade_counts if vv == v)
        if decs:
            row = "  ".join(f"1e{d}:{decade_counts[(v,d)]}" for d in decs)
            print(f"  {v:<14s} {row}")

    if len(lib_rows) > args.max_rows:
        idx = np.random.default_rng(1).choice(
            len(lib_rows), args.max_rows, replace=False)
        lib_rows = [lib_rows[i] for i in sorted(idx)]

    with open(args.outfile, "w") as out:
        out.write("volume,x_mm,y_mm,z_mm\n")
        for v, x, y, z in lib_rows:
            out.write(f"{v},{x:.3f},{y:.3f},{z:.3f}\n")
    print(f"\nLibrary: {len(lib_rows):,} wall-capture vertices → {args.outfile}")

    if args.gas_lib:
        if len(gas_rows) > args.max_rows:
            idx = np.random.default_rng(2).choice(
                len(gas_rows), args.max_rows, replace=False)
            gas_rows = [gas_rows[i] for i in sorted(idx)]
        with open(args.gas_lib, "w") as out:
            out.write("volume,x_mm,y_mm,z_mm\n")
            for v, x, y, z in gas_rows:
                out.write(f"{v},{x:.3f},{y:.3f},{z:.3f}\n")
        print(f"Gas library: {len(gas_rows):,} He3Gas vertices "
              f"(E_n <= {args.gas_emax} eV) → {args.gas_lib}")
    print("Gamma-source event weight = (wall captures/neutron) × (Σ Iγ of the")
    print("cascade table) / N_gamma_generated — see PLAN_NEUTRON_CAMPAIGN.md.")


if __name__ == "__main__":
    main()
