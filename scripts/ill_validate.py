#!/usr/bin/env python3
"""
ill_validate.py — V0/V1/C1 checks on ILL neutron runs (HANDOFF_SIM.md §6)
=========================================================================
Reads EventTree from `mx17_full_sim --beam ill ...` output and prints:

  * generator: accepted / thrown rays, λ moments (particle mean, capture/particle);
  * terminal-interaction budget per (volume, process), per primary neutron;
  * He3Gas absorption: stop-depth quantiles along the beam (cap_y − y_w),
    against the analytic law (x17_facility_search/ill/out/cell_depth.csv);
  * the beam spot at the entrance-window plane (straight-line projection of
    the primary from the gun), against a uniform Ø2a disk;
  * --slab: uncollided transmission (no interaction in "Slab").

Weights: the ³He(n,γ) bias gives captures weight 1/factor; budget rows are
summed with EventTree.weight, i.e. per *unbiased* primary.  (n,p) rows carry
weight 1 at any bias.

    python3 ill_validate.py <files or dir> [--yw -21.6] [--depth-csv cell_depth.csv
        --p-bar 1 --L-cm 30] [--json out.json]
"""

import argparse
import glob
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import uproot


def collect(inputs):
    files = []
    for inp in inputs:
        p = Path(inp)
        if p.is_dir():
            files += sorted(str(f) for f in p.glob("*.root"))
        else:
            files += sorted(glob.glob(inp))
    return files


def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def load(files):
    keys = ["neutron_E_eV", "n_thrown", "capture_vol", "capture_proc", "cap_x", "cap_y",
            "cap_z", "weight", "vtx_x", "vtx_y", "vtx_z", "em_px", "em_py", "em_pz",
            "first_vol", "first_proc", "first_y"]
    parts = defaultdict(list)
    for f in files:
        with uproot.open(f) as fh:
            if "EventTree" not in fh:
                print(f"  skip (no EventTree): {f}")
                continue
            t = fh["EventTree"]
            have = [k for k in keys if k in t.keys()]
            arr = t.arrays(have, library="np")
            for k in have:
                parts[k].append(arr[k])
    out = {k: np.concatenate(v) for k, v in parts.items()}
    for k in ("capture_vol", "capture_proc", "first_vol", "first_proc"):
        if k in out:
            out[k] = dec(out[k])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--yw", type=float, default=None, help="entrance-window y [mm]")
    ap.add_argument("--beam-radius", type=float, default=10.0)
    ap.add_argument("--depth-csv", default=None)
    ap.add_argument("--p-bar", type=float, default=None)
    ap.add_argument("--L-cm", type=float, default=None)
    ap.add_argument("--slab", action="store_true")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    files = collect(a.inputs)
    d = load(files)
    N = len(d["neutron_E_eV"])
    res = dict(files=len(files), n_primaries=N)
    print(f"{len(files)} files, {N:,} primaries")

    # ── generator ───────────────────────────────────────────────────────────
    lam = np.sqrt(81.804e-3 / d["neutron_E_eV"])
    thrown = d["n_thrown"].sum() if "n_thrown" in d else N
    res.update(accept=N / thrown, lam_mean=float(lam.mean()),
               cap_over_part=float((lam / 1.8).mean()),
               lam_p10=float(np.percentile(lam, 10)), frac_below_2A=float((lam < 2).mean()))
    print(f"\nGenerator: accepted/thrown = {N:,}/{thrown:,} = {N/thrown:.4f}")
    print(f"  λ particle mean {lam.mean():.3f} Å (expect 4.9), capture/particle "
          f"{(lam/1.8).mean():.3f} (expect 2.71), 10th pct {np.percentile(lam,10):.2f} Å, "
          f"<2 Å {100*(lam<2).mean():.1f} % (expect 10)")

    # ── spot at the window plane ────────────────────────────────────────────
    if a.yw is not None and "em_py" in d:
        s = (a.yw - d["vtx_y"]) / d["em_py"]
        x = d["vtx_x"] + s * d["em_px"]
        z = d["vtx_z"] + s * d["em_pz"]
        r = np.hypot(x, z)
        R = a.beam_radius
        q = np.percentile(r, [50, 90, 99, 99.9])
        res.update(spot_r50=q[0], spot_r90=q[1], spot_r99=q[2], spot_r999=q[3],
                   spot_out_aperture=float((r > R).mean()),
                   spot_sx=float(x.std()), spot_sz=float(z.std()))
        print(f"\nSpot at y_w = {a.yw} mm: r50/90/99/99.9 = " + "/".join(f"{v:.2f}" for v in q)
              + f" mm (uniform disk r={R}: {R*np.sqrt(.5):.2f}/{R*np.sqrt(.9):.2f}/"
              f"{R*np.sqrt(.99):.2f}); beyond r={R}: {100*(r>R).mean():.2f} %")
        print(f"  σx = {x.std():.2f}, σz = {z.std():.2f} mm (disk: {R/2:.2f})")

    # ── budget ──────────────────────────────────────────────────────────────
    w = d.get("weight", np.ones(N))
    vol, proc = d["capture_vol"], d["capture_proc"]
    vol = np.where(vol == "", "(escaped)", vol)
    print("\nTerminal interaction of the primary, per primary (weighted):")
    budget = {}
    keys = sorted(set(zip(vol, proc)))
    rows = []
    for v, p in keys:
        m = (vol == v) & (proc == p)
        rows.append((w[m].sum() / N, m.sum(), v, p))
    for frac, n, v, p in sorted(rows, reverse=True):
        print(f"  {v:22s} {p:20s} {frac:.4e}  ({n:,} raw)")
        budget[f"{v}|{p}"] = dict(per_n=frac, raw=int(n))
    res["budget"] = budget
    gas_np = (vol == "He3Gas") & (proc == "neutronInelastic")
    res["absorbed_He3_np"] = float(gas_np.mean())
    print(f"  → absorbed in He3Gas by (n,p): {gas_np.mean():.5f}")

    # ── stop depth ──────────────────────────────────────────────────────────
    if a.yw is not None and gas_np.any():
        depth = (d["cap_y"][gas_np] - a.yw) / 10.0     # cm
        qs = [0.10, 0.16, 0.50, 0.84, 0.90, 0.99]
        g = np.quantile(depth, qs)
        res["depth_cm"] = {f"q{int(100*q)}": float(v) for q, v in zip(qs, g)}
        res["depth_mean_cm"] = float(depth.mean())
        rr = np.hypot(d["cap_x"][gas_np], d["cap_z"][gas_np])
        res["vertex_sx_mm"] = float(d["cap_x"][gas_np].std())
        res["vertex_sz_mm"] = float(d["cap_z"][gas_np].std())
        line = "  Geant4 : " + "  ".join(f"q{int(100*q)} {v:6.2f}" for q, v in zip(qs, g)) \
            + f"  mean {depth.mean():.2f} cm"
        print("\nHe3Gas (n,p) stop depth [cm] (from the window):")
        print(line)
        if a.depth_csv and a.p_bar is not None:
            c = np.genfromtxt(a.depth_csv, delimiter=",", names=True)
            m = (np.isclose(c["p_bar"], a.p_bar)) & (np.isclose(c["L_cm"], a.L_cm))
            y, cdf = c["depth_cm"][m], c["cdf"][m]
            an = np.interp(qs, cdf, y)
            print("  analytic: " + "  ".join(f"q{int(100*q)} {v:6.2f}" for q, v in zip(qs, an)))
            res["depth_analytic_cm"] = {f"q{int(100*q)}": float(v) for q, v in zip(qs, an)}
            # KS distance against the analytic CDF
            ds = np.sort(depth)
            ks = np.max(np.abs(np.arange(1, len(ds) + 1) / len(ds) - np.interp(ds, y, cdf)))
            res["depth_ks"] = float(ks)
            print(f"  KS distance {ks:.4f}  (n={len(ds):,}; 95 % crit ≈ {1.36/np.sqrt(len(ds)):.4f})")
        print(f"  vertex σx {res['vertex_sx_mm']:.2f}, σz {res['vertex_sz_mm']:.2f} mm, "
              f"r90 {np.percentile(rr, 90):.2f} mm")

    # ── slab transmission ───────────────────────────────────────────────────
    if a.slab:
        # condition on primaries that reach the slab uncollided: drop those
        # whose first interaction is upstream of it (air between gun and slab)
        fv, fp, fy = d["first_vol"], d["first_proc"], d["first_y"]
        upstream = (fv != "") & (fv != "Slab") & (fy < 0)
        n0 = (~upstream).sum()
        hit = fv == "Slab"
        T = 1 - hit.sum() / n0
        err = np.sqrt(T * (1 - T) / n0)
        res.update(slab_T_uncollided=float(T), slab_T_err=float(err),
                   slab_upstream_scatter=float(upstream.mean()))
        print(f"\nSlab: uncollided transmission {T:.5f} ± {err:.5f}  "
              f"(of {n0:,} reaching it; {100*upstream.mean():.2f} % scattered upstream)")
        for p in sorted(set(fp[hit])):
            print(f"  first interaction in slab: {p:20s} {(fp[hit] == p).sum()/n0:.5f}")

    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
