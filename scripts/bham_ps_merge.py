#!/usr/bin/env python3
"""
bham_ps_merge.py -- merge one wall stage of the Birmingham room runs (submit_room.py --mode ps).

    python3 bham_ps_merge.py <stage dir on EOS> -o <next_gun.csv> [--npz <all_planes.npz>]

Every job wrote <out>.ps.csv (rows: plane, type, E_MeV, x, y, z, dx, dy, dz, w, uncol; w = track × event
weight) and <out>.nev (events).  The stage estimate of any tally per *source particle* is Σ w / N_events,
so each row's weight is divided by the stage's total N.  The rows of the last (killing) plane become the
next stage's --ps-gun file (the generator multiplies by N_rec and draws uniformly, which keeps the
normalisation per source particle).  All planes go to the npz (depth curves, spectra).
"""
import argparse
import glob
import os

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("stage")
ap.add_argument("-o", "--out", required=True)
ap.add_argument("--npz")
a = ap.parse_args()

nev = 0
for f in glob.glob(os.path.join(a.stage, "*.nev")):
    nev += int(open(f).read().split()[0])
parts = []
for f in sorted(glob.glob(os.path.join(a.stage, "*.ps.csv"))):
    if not os.path.exists(f[:-len(".ps.csv")] + ".nev"):
        continue                       # job did not finish: its events are not in nev either
    d = np.loadtxt(f, delimiter=",", skiprows=1, ndmin=2)
    if d.size:
        parts.append(d)
D = np.concatenate(parts) if parts else np.zeros((0, 11))
D[:, 9] /= max(nev, 1)
last = int(D[:, 0].max()) if len(D) else 0
L = D[D[:, 0] == last]
with open(a.out, "w") as f:
    f.write(f"# merged from {a.stage}: N_events {nev}, plane {last}; w per source particle\n")
    f.write("type,E_MeV,x,y,z,dx,dy,dz,w,uncol\n")
    for r in L:
        f.write("%d,%.6g,%.6g,%.6g,%.6g,%.7g,%.7g,%.7g,%.6g,%d\n" % (r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10]))
if a.npz:
    np.savez_compressed(a.npz, plane=D[:, 0].astype(np.int8), type=D[:, 1].astype(np.int8), E=D[:, 2], x=D[:, 3],
                        y=D[:, 4], z=D[:, 5], dx=D[:, 6], dy=D[:, 7], dz=D[:, 8], w=D[:, 9],
                        uncol=D[:, 10].astype(np.int8), n_events=np.array([nev]))
for p in range(last + 1):
    for t, nm in ((0, "n"), (1, "g")):
        m = (D[:, 0] == p) & (D[:, 1] == t)
        print(f"plane {p} {nm}: {m.sum():8d} rows, {D[m, 9].sum():.4g} per source particle "
              f"(uncollided {D[m & (D[:, 10] == 1), 9].sum():.3g})")
print(f"N_events {nev}; {len(L)} rows -> {a.out}")
