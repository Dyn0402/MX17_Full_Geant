#!/usr/bin/env python3
"""Full mechanism characterization of SiPM^plastic legs from the gamma-source
truth run (one Al capture-cascade gamma per event, isotropic, from real
capture vertices; HitTree carries origin_vol/origin_proc/origin_ke/ox,oy,oz).

For every leg (per arm: max single SiPM-bar edep >= 0.5 MIP AND max plastic-bar
edep >= 0.5 MIP) classify:
  - topology: same-track / same-parent / independent
  - for same-track: birth volume (grouped), creator process, birth KE,
    particle species, gamma line energy, birth position
Writes JSON summary + npz histograms.

Usage: python3 analyze_gsrc_mechanism.py <files...> out.json out.npz
"""
import sys, glob, json
import numpy as np
import uproot

files = []
for a in sys.argv[1:-2]:
    files += sorted(glob.glob(a))
OUTJ, OUTN = sys.argv[-2], sys.argv[-1]

MIP_S, MIP_P = 0.458, 4.334 * 2.0 / 2.5
THR_S, THR_P = 0.5 * MIP_S, 0.5 * MIP_P
SIPM, PLAS = ("PlasticScint",), ("BackScintL", "BackScintR")

VOLGROUP = {}
def volgroup(v):
    if v in VOLGROUP: return VOLGROUP[v]
    if v in ("He3Cap_Al", "He3Gas"): g = "Al capsule"
    elif v == "He3Cap_CFRP": g = "capsule CFRP"
    elif v == "World": g = "air"
    elif v.startswith(("GasWindow", "Drift", "Micromesh", "AmpGas",
                       "ResistivePaste", "PCB_")): g = "MM+PCB"
    elif v == "PlasticScint": g = "SiPM bar"
    elif v.startswith("BackScint"): g = "plastic+wrap"
    elif v.startswith(("LiqScint", "LS_")): g = "LS"
    else: g = "other:" + v
    VOLGROUP[v] = g
    return g

def dec(a):
    return np.asarray([x.decode() if isinstance(x, bytes) else x for x in a],
                      dtype=object)

n_ev = 0
n_leg = 0
topo = {"same_track": 0, "same_parent": 0, "independent": 0}
by_volproc = {}          # (volgroup, proc) -> count [same-track legs]
by_species = {}
lineE_leg, lineE_all = [], []
oke, o_r, o_w = [], [], []      # birth KE, dist from beam axis, |pos| along arm
both_pair = 0                    # same-parent e+ AND e- both cross the two dets

for fp in files:
    with uproot.open(fp) as f:
        n_ev += f["EventTree"].num_entries
        et = f["EventTree"].arrays(["eventID", "inv_mass"], library="np")
        lineE = np.zeros(int(et["eventID"].max()) + 1)
        lineE[et["eventID"]] = et["inv_mass"]
        lineE_all.append(et["inv_mass"])
        # pass 1 (chunked, light branches): find legs
        parts = []
        for t in f["HitTree"].iterate(
                ["eventID", "armID", "detType", "edep", "u"],
                library="np", step_size="300 MB"):
            det0 = t["detType"].astype(str)
            k0 = np.isin(det0, SIPM) | np.isin(det0, PLAS)
            if k0.any():
                parts.append({"ev": t["eventID"][k0], "arm": t["armID"][k0],
                              "det": det0[k0], "edep": t["edep"][k0] * 1e-6,
                              "u": t["u"][k0]})
        ev = np.concatenate([p["ev"] for p in parts])
        arm = np.concatenate([p["arm"] for p in parts])
        detname = np.concatenate([p["det"] for p in parts])
        edep0 = np.concatenate([p["edep"] for p in parts])
        u0 = np.concatenate([p["u"] for p in parts])
        sflag = np.isin(detname, SIPM)
        chan = np.where(sflag, np.clip((u0 + 250.) // 25., 0, 19),
                        np.where(detname == "BackScintL", 20, 21)).astype(int)
        key = (ev.astype(np.int64) * 4 + arm) * 22 + chan
        uk, inv = np.unique(key, return_inverse=True)
        bar_edep = np.zeros(len(uk)); np.add.at(bar_edep, inv, edep0)
        ea = (uk // 22)
        uea, inv2 = np.unique(ea, return_inverse=True)
        isbar_s = (uk % 22) < 20
        Smax = np.zeros(len(uea)); Pmax = np.zeros(len(uea))
        np.maximum.at(Smax, inv2[isbar_s], bar_edep[isbar_s])
        np.maximum.at(Pmax, inv2[~isbar_s], bar_edep[~isbar_s])
        leg_ea = uea[(Smax >= THR_S) & (Pmax >= THR_P)]
        leg_events = np.unique(leg_ea // 4)

        # pass 2 (chunked, full branches): keep only hits of leg events
        h_parts = []
        for t in f["HitTree"].iterate(
                ["eventID", "armID", "detType", "trackID", "parentID",
                 "particle", "edep", "u", "origin_vol", "origin_proc",
                 "origin_ke", "ox", "oy", "oz"],
                library="np", step_size="300 MB"):
            m = np.isin(t["eventID"], leg_events)
            if m.any():
                h_parts.append({k: v[m] for k, v in t.items()})
    H = {k: np.concatenate([p[k] for p in h_parts]) for k in h_parts[0]}
    det = dec(H["detType"]); par = dec(H["particle"])
    ovol = dec(H["origin_vol"]); oproc = dec(H["origin_proc"])
    is_s = np.isin(det, SIPM); is_p = np.isin(det, PLAS)
    keep = is_s | is_p
    ev, arm = H["eventID"], H["armID"]
    tid, pid = H["trackID"], H["parentID"]
    edep = H["edep"] * 1e-6            # MeV
    sflag = is_s
    prt, ov, op = par, ovol, oproc
    okev = H["origin_ke"]
    ox, oy, oz = H["ox"], H["oy"], H["oz"]

    ea_hits = np.where(keep, (ev.astype(np.int64) * 4 + arm), -1)
    horder = np.argsort(ea_hits, kind="stable")
    sea = ea_hits[horder]
    for x in leg_ea:
        n_leg += 1
        lo = np.searchsorted(sea, x, "left")
        hi = np.searchsorted(sea, x, "right")
        idx = horder[lo:hi]
        sf = sflag[idx]
        t_i, p_i = tid[idx], pid[idx]
        ts, tp = set(t_i[sf].tolist()), set(t_i[~sf].tolist())
        ps, pp = set(p_i[sf].tolist()), set(p_i[~sf].tolist())
        common = ts & tp
        e = int(x // 4)
        lineE_leg.append(float(lineE[e]))
        if common:
            topo["same_track"] += 1
            # pick the common track with the largest plastic deposit
            best, bde = None, -1.0
            for t in common:
                de = edep[idx[(~sf) & (t_i == t)]].sum()
                if de > bde: bde, best = de, t
            j = idx[t_i == best][0]
            k = (volgroup(str(ov[j])), str(op[j]))
            by_volproc[k] = by_volproc.get(k, 0) + 1
            sp = str(prt[j]); by_species[sp] = by_species.get(sp, 0) + 1
            oke.append(float(okev[j]))
            o_r.append(float(np.hypot(ox[j], oz[j])))
            o_w.append(float(oy[j]))
            # pair signature: e+ and e- siblings both crossing
            if len(common) >= 2:
                sps = {str(prt[idx[t_i == t][0]]) for t in common}
                if "e+" in sps and "e-" in sps:
                    both_pair += 1
        elif (ps & pp) or (ts & pp) or (ps & tp):
            topo["same_parent"] += 1
        else:
            topo["independent"] += 1
    print(f"done {fp.split('/')[-1]}: cum legs={n_leg}", flush=True)

lineE_all = np.concatenate(lineE_all)
out = {
    "n_events": int(n_ev), "n_legs": int(n_leg), "topology": topo,
    "legs_per_gamma": n_leg / n_ev,
    "same_track_by_volproc": {f"{v}|{p}": c for (v, p), c in
                              sorted(by_volproc.items(), key=lambda kv: -kv[1])},
    "species": by_species, "both_pair_cross": both_pair,
    "origin_ke_median": float(np.median(oke)) if oke else None,
}
json.dump(out, open(OUTJ, "w"), indent=1)
np.savez(OUTN, lineE_leg=np.array(lineE_leg), lineE_all=lineE_all,
         oke=np.array(oke), o_r=np.array(o_r), o_w=np.array(o_w))
print(json.dumps(out, indent=1))
