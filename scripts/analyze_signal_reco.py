#!/usr/bin/env python3
"""analyze_signal_reco.py — Task 3/4c for the SIGNAL: MM 2-track acceptance and
reconstructed (MSC-smeared) opening angle for X17 and IPC, from pairs_thermal_
trig_2cm, using the SAME MM-hit reconstruction as the Al background
(analyze_convpair_reco.py). This makes the Al-vs-IPC-vs-X17 comparison
apples-to-apples: reconstructed opening angle among MM 2-track pairs.

Signal leptons are the two PRIMARIES (parentID==0), split e-/e+ by `particle`.
MM-entry direction = momentum dir of earliest-time DriftGas hit.

Writes /afs/cern.ch/work/d/dneff/convpair_out/signal_reco.npz.
Run on lxplus.  Optional arg NFILES (default all).
"""
import sys, glob
import numpy as np, uproot

PAIRS = "/eos/experiment/ntof/data/x17/full_sim/pairs_thermal_trig_2cm"
OUT   = "/afs/cern.ch/work/d/dneff/convpair_out/signal_reco.npz"
NF = int(sys.argv[1]) if len(sys.argv) > 1 else 100

def process(fp, rows):
    f = uproot.open(fp)
    e = f["EventTree"].arrays(
        ["eventID", "event_type", "openingAngle", "em_ke", "ep_ke"], library="np")
    et = {int(ev): (int(t), oa, mn) for ev, t, oa, mn in zip(
        e["eventID"], e["event_type"], e["openingAngle"],
        np.minimum(e["em_ke"], e["ep_ke"]))}
    # accumulate per-chunk earliest-hit candidates per (event,lepton) key.
    # key = eventID*2 + (0:e-, 1:e+). Vectorized reduction — no per-hit loop.
    cand = []   # list of (key, time, px,py,pz, arm)
    for t in f["HitTree"].iterate(
            ["eventID", "parentID", "particle", "detType", "armID",
             "px", "py", "pz", "time"], library="np", step_size="300 MB"):
        prim = t["parentID"] == 0
        if not prim.any():
            continue
        det = t["detType"]
        gasi = prim & ((det == b"DriftGas") | (det == "DriftGas"))
        if not gasi.any():
            continue
        idx = np.where(gasi)[0]
        part = t["particle"][idx]
        code = np.where((part == b"e-") | (part == "e-"), 0,
                        np.where((part == b"e+") | (part == "e+"), 1, -1))
        good = code >= 0
        idx = idx[good]; code = code[good]
        if len(idx) == 0:
            continue
        key = t["eventID"][idx].astype(np.int64) * 2 + code
        tm = t["time"][idx]
        order = np.lexsort((tm, key))
        ks = key[order]; first = np.concatenate(([True], ks[1:] != ks[:-1]))
        fi = order[first]
        cand.append(np.column_stack([key[fi], tm[fi], t["px"][idx][fi],
                                     t["py"][idx][fi], t["pz"][idx][fi],
                                     t["armID"][idx][fi].astype(float)]))
    best = {}  # key -> (time, px,py,pz, arm)
    if cand:
        C = np.concatenate(cand)
        order = np.lexsort((C[:, 1], C[:, 0]))
        Cs = C[order]; k = Cs[:, 0]
        firstmask = np.concatenate(([True], k[1:] != k[:-1]))
        for row in Cs[firstmask]:
            best[int(row[0])] = row[1:]
    for ev, (typ, oa, emn) in et.items():
        rec = {"type": typ, "theta_truth": oa, "emin_truth": emn,
               "theta_reco": np.nan, "n_mm": 0, "same_arm": -1}
        dirs = {}; arms = {}
        for code, lep in ((0, "e-"), (1, "e+")):
            r = best.get(ev * 2 + code)
            if r is None:
                continue
            v = r[1:4]; n = np.linalg.norm(v)
            if n == 0:
                continue
            dirs[lep] = v / n; arms[lep] = int(r[4])
        rec["n_mm"] = len(dirs)
        if len(dirs) == 2:
            coso = float(np.clip(np.dot(dirs["e-"], dirs["e+"]), -1, 1))
            rec["theta_reco"] = np.degrees(np.arccos(coso))
            rec["same_arm"] = int(arms["e-"] == arms["e+"])
        rows.append(rec)

def main():
    files = sorted(glob.glob(f"{PAIRS}/*.root"))[:NF]
    rows = []
    for k, fp in enumerate(files):
        process(fp, rows)
        print(f"  [{k+1}/{len(files)}] {fp.split('/')[-1]}  cum={len(rows)}", flush=True)
    typ = np.array([r["type"] for r in rows])
    tt  = np.array([r["theta_truth"] for r in rows])
    tr  = np.array([r["theta_reco"] for r in rows])
    nmm = np.array([r["n_mm"] for r in rows])
    emn = np.array([r["emin_truth"] for r in rows])
    same = np.array([r["same_arm"] for r in rows])
    for name, code in [("X17", 0), ("IPC", 1)]:
        m = typ == code; pair = m & (nmm == 2)
        print(f"\n{name}: generated={m.sum()}  MM 2-track acceptance={pair.sum()/m.sum()*100:.2f}%")
        if pair.any():
            print(f"  truth θ (MM pairs): median={np.median(tt[pair]):.1f}  >108°={(tt[pair]>108).mean()*100:.1f}%")
            print(f"  reco  θ (MM pairs): median={np.median(tr[pair]):.1f}  >108°={(tr[pair]>108).mean()*100:.1f}%")
            print(f"  same-arm: {(same[pair]==1).mean()*100:.1f}%")
    np.savez(OUT, type=typ, theta_truth=tt, theta_reco=tr, n_mm=nmm,
             emin_truth=emn, same_arm=same)
    print(f"\nwrote {OUT}")

if __name__ == "__main__":
    main()
