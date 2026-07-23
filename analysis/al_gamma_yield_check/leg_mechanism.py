#!/usr/bin/env python3
"""Classify what makes a per-arm SiPM^plastic coincidence ('leg', 0.5 MIP each)
in the thermal campaign. For each leg: do the SiPM-wall hits and the
BackScint hits share a track (charged particle crossing both), share a parent
(same gamma), or come from independent particles (e.g. different cascade
gammas)?  Also tallies capture_vol of the event and particle species.

Run on lxplus from the repo root:
  python3 leg_mechanism.py /eos/.../neutrons_thermal_trig_2cm NF out.json
"""
import sys, glob, json
import numpy as np
import uproot

sys.path.insert(0, "scripts")
from analyze_trigger_thermal import digest_file

DIR = sys.argv[1]
NF = int(sys.argv[2]) if len(sys.argv) > 2 else 4
OUT = sys.argv[3] if len(sys.argv) > 3 else "leg_mechanism.json"

MIP_S, MIP_P = 0.458, 4.334 * 2.0 / 2.5
FILES = sorted(glob.glob(f"{DIR}/*.root"))[:NF]

SIPM, PLAS = ("PlasticScint",), ("BackScintL", "BackScintR")

cnt = {"legs": 0, "events": 0, "nev_tot": 0,
       "same_track": 0, "same_parent": 0, "independent": 0,
       "st_with_driftgas": 0,      # crossing track also crossed the MM drift
       "st_pair_both_cross": 0,    # e+ AND e- siblings both cross both dets
       "cap_He3Cap_Al": 0, "cap_other": {}}
species = {}          # particle of the SiPM-hit track for same_track legs
ke_same = []          # ke [MeV] of that track at the SiPM hit

for fp in FILES:
    obs, aux, nev = digest_file(fp, gate=True)
    cnt["nev_tot"] += nev
    S = obs[:, :, 0] / MIP_S
    P = obs[:, :, 1] / MIP_P
    leg = (S >= 0.5) & (P >= 0.5)                      # (N,4)
    sel = leg.any(axis=1)
    ev_leg = aux["eventID"][sel]
    arms_leg = leg[sel]
    if not len(ev_leg):
        continue
    ev_set = set(int(e) for e in ev_leg)

    with uproot.open(fp) as f:
        et = f["EventTree"].arrays(["eventID", "capture_vol"], library="np")
        cap = dict(zip(et["eventID"].tolist(), et["capture_vol"].tolist()))
        h_parts = []
        for t in f["HitTree"].iterate(
                ["eventID", "armID", "detType", "trackID", "parentID",
                 "particle", "edep", "ke"], library="np", step_size="200 MB"):
            m = np.isin(t["eventID"], ev_leg)
            if m.any():
                h_parts.append({k: v[m] for k, v in t.items()})
    h = {k: np.concatenate([p[k] for p in h_parts]) for k in h_parts[0]}

    for ev, arms in zip(ev_leg, arms_leg):
        cnt["events"] += 1
        cv = cap.get(int(ev), "")
        if cv == "He3Cap_Al":
            cnt["cap_He3Cap_Al"] += 1
        else:
            cnt["cap_other"][cv] = cnt["cap_other"].get(cv, 0) + 1
        em = h["eventID"] == ev
        for a in np.nonzero(arms)[0]:
            cnt["legs"] += 1
            am = em & (h["armID"] == a)
            is_s = am & np.isin(h["detType"], SIPM)
            is_p = am & np.isin(h["detType"], PLAS)
            ts, tp = set(h["trackID"][is_s].tolist()), set(h["trackID"][is_p].tolist())
            ps, pp = set(h["parentID"][is_s].tolist()), set(h["parentID"][is_p].tolist())
            common = ts & tp
            if common:
                cnt["same_track"] += 1
                tid = next(iter(common))
                i = np.nonzero(is_s & (h["trackID"] == tid))[0]
                if len(i):
                    sp = str(h["particle"][i[0]])
                    species[sp] = species.get(sp, 0) + 1
                    ke_same.append(float(h["ke"][i[0]]))
                if ((em & (h["detType"] == "DriftGas")
                        & (h["trackID"] == tid)).any()):
                    cnt["st_with_driftgas"] += 1
                if len(common) >= 2:
                    sps = set()
                    for t in common:
                        j = np.nonzero(em & (h["trackID"] == t))[0]
                        if len(j):
                            sps.add(str(h["particle"][j[0]]))
                    if "e+" in sps and "e-" in sps:
                        cnt["st_pair_both_cross"] += 1
            elif (ps & pp) or (ts & pp) or (ps & tp):
                cnt["same_parent"] += 1
            else:
                cnt["independent"] += 1
    print(f"done {fp.split('/')[-1]}  cum legs={cnt['legs']}", flush=True)

cnt["species_same_track"] = species
cnt["ke_same_track_median"] = float(np.median(ke_same)) if ke_same else None
cnt["ke_same_track_q"] = ([float(q) for q in np.percentile(ke_same, [25, 75])]
                          if ke_same else None)
json.dump(cnt, open(OUT, "w"), indent=1)
print(json.dumps(cnt, indent=1))
