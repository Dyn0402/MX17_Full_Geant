#!/usr/bin/env python3
"""count_singles_thermal.py — in-gate singles multiplicity per pulse.
Counts individual detector bars firing above 0.5 MIP (the 'coincident singles'
seen in a readout), from raw thermal neutron files. SiPM MIP=458 keV (geometry
unchanged, robust); plastic MIP quoted for both 2.0 and 2.5 cm to bracket."""
import glob, sys
import numpy as np, uproot

DIR = sys.argv[1] if len(sys.argv) > 1 else \
      "/eos/experiment/ntof/data/x17/full_sim/neutrons_thermal_trig_2cm"
NF  = int(sys.argv[2]) if len(sys.argv) > 2 else 3
FILES = sorted(glob.glob(f"{DIR}/*.root"))[:NF]

MIP_SIPM = 0.458            # MeV/MIP (3 mm SiPM bar, unchanged)
MIP_PL25, MIP_PL20 = 4.334, 4.334 * 2.0 / 2.5   # 2.5 cm / scaled 2.0 cm
FLIGHT, GATE, N_PULSE = 1.41, 1.0, 4.284e6

def codes(det):
    d = np.asarray(det).astype(str); c = np.full(len(d), -1, np.int8)
    c[np.char.startswith(d, "PlasticScint")] = 0
    c[np.char.startswith(d, "BackScintL")] = 1
    c[np.char.startswith(d, "BackScintR")] = 2
    return c

tS = tP20 = tP25 = 0; nev_tot = 0
for fp in FILES:
    with uproot.open(fp) as f:
        nev = f["EventTree"].num_entries; nev_tot += nev
        et = f["EventTree"].arrays(["eventID", "neutron_E_eV"], library="np")
        eN = np.zeros(int(et["eventID"].max()) + 1); eN[et["eventID"]] = et["neutron_E_eV"]
        for t in f["HitTree"].iterate(["eventID","armID","detType","edep","time","u"],
                                       library="np", step_size="300 MB"):
            c = codes(t["detType"]); keep = c >= 0
            e_hit = eN[t["eventID"]]
            with np.errstate(divide="ignore"):
                tof = np.where(e_hit > 0, FLIGHT / np.sqrt(e_hit), np.inf)
            keep &= (tof + t["time"] * 1e-6) > GATE          # in-gate (>1 ms)
            c = c[keep]; ev = t["eventID"][keep].astype(np.int64)
            arm = t["armID"][keep].astype(np.int64); edep = t["edep"][keep] * 1e-6
            u = t["u"][keep]
            chan = np.full(len(ev), 22, np.int64)
            s = c == 0; chan[s] = np.clip((u[s] + 250.) // 25., 0, 19).astype(np.int64)
            chan[c == 1] = 20; chan[c == 2] = 21
            key = (ev * 4 + arm) * 23 + chan
            uk, inv = np.unique(key, return_inverse=True)
            es = np.zeros(len(uk)); np.add.at(es, inv, edep)
            kch = uk % 23
            tS   += int(np.sum((kch <= 19) & (es > 0.5 * MIP_SIPM)))
            tP20 += int(np.sum(((kch == 20) | (kch == 21)) & (es > 0.5 * MIP_PL20)))
            tP25 += int(np.sum(((kch == 20) | (kch == 21)) & (es > 0.5 * MIP_PL25)))

print(f"files={len(FILES)}  neutrons={nev_tot:,}  N_pulse={N_PULSE:.3e}")
print(f"SiPM bar-singles >0.5 MIP (229 keV) : {tS/nev_tot:.3e}/n -> {tS/nev_tot*N_PULSE:,.0f}/pulse")
print(f"plastic singles  >0.5 MIP @2.0cm(1.73MeV): {tP20/nev_tot*N_PULSE:,.0f}/pulse")
print(f"plastic singles  >0.5 MIP @2.5cm(2.17MeV): {tP25/nev_tot*N_PULSE:,.0f}/pulse")
