#!/usr/bin/env python3
"""extract_signal_openingangle.py — Task 2: IPC & X17 truth opening-angle
(and individual-lepton-energy) distributions, straight from the pairs tree.
Also pulls the low-stats Al conv-pair opening angle from the trajdump as a
cross-check (Task 1a). Writes a small npz for local plotting.

Run on lxplus (needs uproot + EOS access):
    python3 scripts/extract_signal_openingangle.py [NFILES]
"""
import sys, glob, csv
from collections import defaultdict
import numpy as np, uproot

PAIRS = "/eos/experiment/ntof/data/x17/full_sim/pairs_thermal_trig_2cm"
TRAJ  = "/tmp/dneff/trajbig_traj_t0.csv"
NF = int(sys.argv[1]) if len(sys.argv) > 1 else 20
OUT = "/tmp/dneff/signal_openingangle.npz"

# ── Task 2: IPC & X17 truth from the pairs EventTree ────────────────────────
br = ["event_type", "openingAngle", "em_ke", "ep_ke", "inv_mass",
      "vtx_x", "vtx_y", "vtx_z", "neutron_E_eV"]
acc = {k: [] for k in br}
files = sorted(glob.glob(f"{PAIRS}/*.root"))[:NF]
for fp in files:
    with uproot.open(fp) as f:
        a = f["EventTree"].arrays(br, library="np")
        for k in br:
            acc[k].append(a[k])
for k in br:
    acc[k] = np.concatenate(acc[k])
et = acc["event_type"]
isX17 = et == 0
isIPC = et == 1
print(f"pairs files={len(files)}  X17={isX17.sum()}  IPC={isIPC.sum()}")
print(f"  X17 openingAngle median={np.median(acc['openingAngle'][isX17]):.1f} deg")
print(f"  IPC openingAngle median={np.median(acc['openingAngle'][isIPC]):.1f} deg")

# ── Task 1a: Al conv-pair opening angle from trajdump (low-stats cross-check) ─
# Reconstruct e+/e- birth directions from the first step of each daughter of a
# gamma that underwent conv. Direction ~= (post-pre) of the daughter's 1st step.
al_theta = []; al_egamma = []; al_r = []; al_emin_ke = []
try:
    gamma_tracks = defaultdict(set)      # ev -> gamma trackIDs
    first = defaultdict(dict)            # (ev,pid) -> {particle: row(firststep)}
    with open(TRAJ) as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            if len(row) < 16: continue
            ev = int(row[0]); tid = int(row[1]); pid = int(row[2]); part = row[3]
            istep = int(row[4]); vol = row[13]; proc = row[14]
            if part == "gamma":
                gamma_tracks[ev].add(tid)
            elif part in ("e+", "e-"):
                d = first[(ev, pid)]
                if part not in d or istep < int(d[part][4]):
                    d[part] = row
    for (ev, pid), d in first.items():
        if pid not in gamma_tracks[ev] or "e+" not in d or "e-" not in d:
            continue
        rows = {}
        for p in ("e+", "e-"):
            x0, y0, z0 = float(d[p][5]), float(d[p][6]), float(d[p][7])
            x1, y1, z1 = float(d[p][8]), float(d[p][9]), float(d[p][10])
            dvec = np.array([x1 - x0, y1 - y0, z1 - z0])
            n = np.linalg.norm(dvec)
            if n == 0: break
            rows[p] = (dvec / n, float(d[p][11]), (x0, y0, z0), d[p][13])
        if len(rows) != 2: continue
        # restrict to Al capsule conversions
        if rows["e+"][3] != "He3Cap_Al": continue
        coso = np.clip(np.dot(rows["e+"][0], rows["e-"][0]), -1, 1)
        al_theta.append(np.degrees(np.arccos(coso)))
        al_egamma.append(rows["e+"][1] + rows["e-"][1] + 1.022)
        x0, y0, z0 = rows["e+"][2]
        al_r.append(np.hypot(x0, z0))               # transverse radius (beam=+Y)
        al_emin_ke.append(min(rows["e+"][1], rows["e-"][1]))
    print(f"trajdump Al conv pairs: {len(al_theta)}  "
          f"median theta={np.median(al_theta) if al_theta else float('nan'):.1f} deg")
except FileNotFoundError:
    print(f"WARNING: trajdump {TRAJ} gone; skipping Al cross-check")

np.savez(OUT,
         x17_theta=acc["openingAngle"][isX17], ipc_theta=acc["openingAngle"][isIPC],
         x17_em=acc["em_ke"][isX17], x17_ep=acc["ep_ke"][isX17],
         ipc_em=acc["em_ke"][isIPC], ipc_ep=acc["ep_ke"][isIPC],
         x17_mee=acc["inv_mass"][isX17], ipc_mee=acc["inv_mass"][isIPC],
         al_theta=np.array(al_theta), al_egamma=np.array(al_egamma),
         al_r=np.array(al_r), al_emin_ke=np.array(al_emin_ke))
print(f"wrote {OUT}")
