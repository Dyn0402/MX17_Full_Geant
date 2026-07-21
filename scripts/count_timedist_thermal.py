#!/usr/bin/env python3
"""count_timedist_thermal.py — time distributions (arrival time t>1 ms) of
SiPM-wall singles, plastic singles, and per-arm SiPM^plastic coincidences
('legs'), per pulse, from the thermal _2cm data.  Arm-level (= DAQ wall-level:
SiPM full wall OR, plastic OR) at 0.5 MIP.  Writes histograms to an npz for
local plotting.  Time = TOF from birth energy (1.41/sqrt(E_eV) ms, 19.5 m);
in-target moderation adds <~0.5 ms smearing at the cold end.
"""
import glob, sys
import numpy as np, uproot
sys.path.insert(0, "scripts")
from analyze_trigger_thermal import digest_file

DIR = sys.argv[1] if len(sys.argv) > 1 else \
      "/eos/experiment/ntof/data/x17/full_sim/neutrons_thermal_trig_2cm"
NF  = int(sys.argv[2]) if len(sys.argv) > 2 else 8
OUT = sys.argv[3] if len(sys.argv) > 3 else "analysis/thermal_2cm/timedist_2cm.npz"
FILES = sorted(glob.glob(f"{DIR}/*.root"))[:NF]

MIP_S = 0.458
MIP_P20, MIP_P25 = 4.334 * 2.0 / 2.5, 4.334
FL, A, B = 1.41, 0.5, 0.5
tedges = np.arange(1.0, 40.01, 0.5)
tc = 0.5 * (tedges[:-1] + tedges[1:])
keys = ["sipm", "plas20", "plas25", "leg20", "leg25"]
H = {k: np.zeros(len(tc)) for k in keys}
nev = 0
for fp in FILES:
    obs, aux, n = digest_file(fp, gate=True)      # in-gate (>1 ms)
    nev += n
    with uproot.open(fp) as f:
        et = f["EventTree"].arrays(["eventID", "neutron_E_eV"], library="np")
    eN = np.zeros(int(et["eventID"].max()) + 1); eN[et["eventID"]] = et["neutron_E_eV"]
    E = eN[aux["eventID"]]
    with np.errstate(divide="ignore"):
        tof = np.where(E > 0, FL / np.sqrt(E), np.inf)
    S = obs[:, :, 0] / MIP_S
    P20 = obs[:, :, 1] / MIP_P20; P25 = obs[:, :, 1] / MIP_P25
    sa = S >= A; pa20 = P20 >= B; pa25 = P25 >= B
    H["sipm"]   += np.histogram(tof, tedges, weights=sa.sum(1))[0]
    H["plas20"] += np.histogram(tof, tedges, weights=pa20.sum(1))[0]
    H["plas25"] += np.histogram(tof, tedges, weights=pa25.sum(1))[0]
    H["leg20"]  += np.histogram(tof, tedges, weights=(sa & pa20).sum(1))[0]
    H["leg25"]  += np.histogram(tof, tedges, weights=(sa & pa25).sum(1))[0]
    print(f"  {fp.split('/')[-1]}: cum events={nev:,}", flush=True)

np.savez(OUT, tedges=tedges, tc=tc, n_events=nev, n_pulse=4.284e6,
         mip_sipm=MIP_S, mip_plas20=MIP_P20, mip_plas25=MIP_P25, **H)
print(f"wrote {OUT}  nev={nev:,}  "
      f"totals/pulse: sipm={H['sipm'].sum()/nev*4.284e6:.0f} "
      f"plas20={H['plas20'].sum()/nev*4.284e6:.0f} leg20={H['leg20'].sum()/nev*4.284e6:.0f}")
