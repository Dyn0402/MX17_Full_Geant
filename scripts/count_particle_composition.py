#!/usr/bin/env python3
"""count_particle_composition.py — what particle deposits the background energy.
Tallies in-gate energy deposition (and hit counts) in the SiPM-wall and plastic
bars by the depositing particle type, from raw thermal data. Shows whether the
background is gamma-induced (Compton/photo electrons) or neutron-recoil (protons)
etc.  detType classes: 0=SiPM bar, 1/2=plastic L/R.
"""
import glob, sys
import numpy as np, uproot
sys.path.insert(0, "scripts")
from analyze_trigger_thermal import det_class_codes

DIR = sys.argv[1] if len(sys.argv) > 1 else \
      "/eos/experiment/ntof/data/x17/full_sim/neutrons_thermal_trig_2cm"
NF  = int(sys.argv[2]) if len(sys.argv) > 2 else 3
FILES = sorted(glob.glob(f"{DIR}/*.root"))[:NF]
FL, GATE = 1.41, 1.0

PARTS = ["e-", "gamma", "proton", "alpha", "triton", "deuteron", "e+", "neutron"]
PIDX = {p: i for i, p in enumerate(PARTS)}
NOTH = len(PARTS)            # index for "other"
DETS = ["SiPM", "plastic"]

edep = np.zeros((2, NOTH + 1))     # [det, particle] MeV
hits = np.zeros((2, NOTH + 1), dtype=np.int64)
nev = 0
for fp in FILES:
    with uproot.open(fp) as f:
        nev += f["EventTree"].num_entries
        et = f["EventTree"].arrays(["eventID", "neutron_E_eV"], library="np")
        eN = np.zeros(int(et["eventID"].max()) + 1); eN[et["eventID"]] = et["neutron_E_eV"]
        for t in f["HitTree"].iterate(["eventID", "detType", "particle", "edep", "time"],
                                       library="np", step_size="300 MB"):
            code = det_class_codes(t["detType"])
            det = np.where(code == 0, 0, np.where((code == 1) | (code == 2), 1, -1))
            keep = det >= 0
            e_hit = eN[t["eventID"]]
            with np.errstate(divide="ignore"):
                tof = np.where(e_hit > 0, FL / np.sqrt(e_hit), np.inf)
            keep &= (tof + t["time"] * 1e-6) > GATE
            det = det[keep]
            ed  = t["edep"][keep] * 1e-6                    # MeV
            praw = np.asarray([x.decode() if isinstance(x, bytes) else x
                               for x in t["particle"][keep]], dtype=object)
            pcat = np.array([PIDX.get(p, NOTH) for p in praw], dtype=np.int64)
            np.add.at(edep, (det, pcat), ed)
            np.add.at(hits, (det, pcat), 1)

labels = PARTS + ["other"]
print(f"files={len(FILES)} events={nev:,}  (in-gate SiPM+plastic deposits)")
for di, dn in enumerate(DETS):
    tot = edep[di].sum()
    print(f"\n{dn}: total {tot:.3e} MeV in-gate")
    order = np.argsort(-edep[di])
    for j in order:
        if edep[di][j] <= 0:
            continue
        print(f"  {labels[j]:9s}  edep {edep[di][j]/tot*100:5.1f}%   hits {hits[di][j]/hits[di].sum()*100:5.1f}%")
np.savez("analysis/thermal_2cm/particle_composition_2cm.npz",
         edep=edep, hits=hits, labels=np.array(labels), dets=np.array(DETS), n_events=nev)
print("\nwrote analysis/thermal_2cm/particle_composition_2cm.npz")
