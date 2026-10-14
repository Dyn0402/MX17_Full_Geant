#!/usr/bin/env python3
"""
ill_make_lib.py — vertex library (volume,x_mm,y_mm,z_mm) from ILL neutron runs,
for `mx17_full_sim --pair-vertex-lib` (S1: He3Gas (n,p); S2: window captures).

    python3 ill_make_lib.py /eos/.../ill/C1/G1 -o lib_G1_gas.csv
    python3 ill_make_lib.py /eos/.../ill/C1w/G1 -o lib_G1_window.csv \\
        --volume He3Cell_Window --proc nCapture
    python3 ill_make_lib.py ... --shift-y 14.4      # S1p: move the cell by Δy

(n,p) and (n,γ) in ³He both go as 1/v, so the (n,p) positions are the (n,γ)
vertex distribution.  Rows are thinned uniformly to --max-rows.  Biased runs
are fine as long as all selected rows share one weight (true for one volume).
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import uproot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--volume", default="He3Gas")
    ap.add_argument("--proc", default="neutronInelastic")
    ap.add_argument("--max-rows", type=int, default=2_000_000)
    ap.add_argument("--shift-y", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    files = []
    for i in a.inputs:
        p = Path(i)
        files += sorted(str(f) for f in p.glob("*_t0.root")) if p.is_dir() else sorted(glob.glob(i))
    rows, wts = [], []
    for f in files:
        t = uproot.open(f)["EventTree"].arrays(["capture_vol", "capture_proc", "cap_x", "cap_y",
                                                "cap_z", "weight"], library="np")
        v = np.char.rstrip(t["capture_vol"].astype(str), "\x00")
        p = np.char.rstrip(t["capture_proc"].astype(str), "\x00")
        p = np.char.replace(np.char.replace(p, "biasWrapper(", ""), ")", "")
        m = (v == a.volume) & (p == a.proc)
        rows.append(np.stack([t["cap_x"][m], t["cap_y"][m] + a.shift_y, t["cap_z"][m]], 1))
        wts.append(t["weight"][m])
    X, W = np.concatenate(rows), np.concatenate(wts)
    rng = np.random.default_rng(a.seed)
    if len(X) > a.max_rows:
        X = X[rng.choice(len(X), a.max_rows, replace=False)]
    with open(a.out, "w") as fh:
        fh.write(f"# {len(files)} files, volume {a.volume}, proc {a.proc}, shift_y {a.shift_y} mm; "
                 f"weight min/max {W.min():.4g}/{W.max():.4g}\n")
        fh.write("volume,x_mm,y_mm,z_mm\n")
        for x, y, z in X:
            fh.write(f"{a.volume},{x:.3f},{y:.3f},{z:.3f}\n")
    print(f"{a.out}: {len(X):,} rows from {len(files)} files; y median {np.median(X[:, 1]):.2f} mm, "
          f"weights {W.min():.4g}–{W.max():.4g}")


if __name__ == "__main__":
    main()
