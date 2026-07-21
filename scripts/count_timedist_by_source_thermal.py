#!/usr/bin/env python3
"""count_timedist_by_source_thermal.py — like count_timedist_thermal.py but the
SiPM/plastic/leg rates are split by the neutron-capture SOURCE (from capture_vol
of each event's primary neutron; every hit in an event descends from that one
capture). Writes per-source time histograms to npz for a stacked local plot.
"""
import glob, sys
import numpy as np, uproot
sys.path.insert(0, "scripts")
from analyze_trigger_thermal import digest_file

DIR = sys.argv[1] if len(sys.argv) > 1 else \
      "/eos/experiment/ntof/data/x17/full_sim/neutrons_thermal_trig_2cm"
NF  = int(sys.argv[2]) if len(sys.argv) > 2 else 8
OUT = sys.argv[3] if len(sys.argv) > 3 else "analysis/thermal_2cm/timedist_bysource_2cm.npz"
FILES = sorted(glob.glob(f"{DIR}/*.root"))[:NF]

MIP_S, MIP_P = 0.458, 4.334 * 2.0 / 2.5      # 2.0 cm plastic
FL, A, B = 1.41, 0.5, 0.5
tedges = np.arange(1.0, 40.01, 0.5); tc = 0.5 * (tedges[:-1] + tedges[1:])
SRC = ["Al capsule", "Capsule CFRP", "LS + vessel", "Back plastic",
       "SiPM scint", "Gas/other"]
NSRC = len(SRC)

def srccode(vraw):
    v = np.asarray([x.decode() if isinstance(x, bytes) else x for x in vraw],
                   dtype=object).astype(str)
    c = np.full(len(v), 5, dtype=np.int8)           # 5 = gas/other/escaped
    c[np.char.startswith(v, "He3Cap_Al")]   = 0
    c[np.char.startswith(v, "He3Cap_CFRP")] = 1
    c[np.char.startswith(v, "LiqScint")]    = 2
    c[np.char.startswith(v, "LS_")]         = 2
    c[np.char.startswith(v, "BackScint")]   = 3
    c[np.char.startswith(v, "PlasticScint")] = 4
    return c

H = {k: np.zeros((NSRC, len(tc))) for k in ["sipm", "plas", "leg"]}
nev = 0
for fp in FILES:
    obs, aux, n = digest_file(fp, gate=True); nev += n
    with uproot.open(fp) as f:
        et = f["EventTree"].arrays(["eventID", "neutron_E_eV", "capture_vol"],
                                    library="np")
    mid = int(et["eventID"].max())
    eN = np.zeros(mid + 1); eN[et["eventID"]] = et["neutron_E_eV"]
    eS = np.full(mid + 1, 5, np.int8); eS[et["eventID"]] = srccode(et["capture_vol"])
    E = eN[aux["eventID"]]; src = eS[aux["eventID"]]
    with np.errstate(divide="ignore"):
        tof = np.where(E > 0, FL / np.sqrt(E), np.inf)
    S = obs[:, :, 0] / MIP_S; P = obs[:, :, 1] / MIP_P
    sa = S >= A; pa = P >= B
    ns, nplast, nl = sa.sum(1), pa.sum(1), (sa & pa).sum(1)
    for s in range(NSRC):
        m = src == s
        if not m.any():
            continue
        H["sipm"][s] += np.histogram(tof[m], tedges, weights=ns[m])[0]
        H["plas"][s] += np.histogram(tof[m], tedges, weights=nplast[m])[0]
        H["leg"][s]  += np.histogram(tof[m], tedges, weights=nl[m])[0]
    print(f"  {fp.split('/')[-1]}: cum events={nev:,}", flush=True)

np.savez(OUT, tedges=tedges, tc=tc, n_events=nev, n_pulse=4.284e6,
         sources=np.array(SRC), mip_sipm=MIP_S, mip_plas=MIP_P, **H)
w = 4.284e6 / nev
print(f"wrote {OUT}  nev={nev:,}")
for s in range(NSRC):
    print(f"  {SRC[s]:16s} legs/pulse={H['leg'][s].sum()*w:7.1f}  "
          f"sipm/pulse={H['sipm'][s].sum()*w:7.0f}  plastic/pulse={H['plas'][s].sum()*w:7.0f}")
