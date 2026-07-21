#!/usr/bin/env python3
"""count_legs_thermal.py — per-arm SiPM AND plastic coincidences ('legs') per
pulse, in-gate. A 'leg' = an arm with (max SiPM bar >0.5 MIP) AND (max plastic
bar >0.5 MIP) — i.e. what the DAQ sees as a single-arm SiPM+plastic single.
Uses the trigger analysis's own per-arm digest so the definition is identical.
"""
import glob, sys
import numpy as np
sys.path.insert(0, "scripts")
from analyze_trigger_thermal import digest_file

DIR = sys.argv[1] if len(sys.argv) > 1 else \
      "/eos/experiment/ntof/data/x17/full_sim/neutrons_thermal_trig_2cm"
NF  = int(sys.argv[2]) if len(sys.argv) > 2 else 3
FILES = sorted(glob.glob(f"{DIR}/*.root"))[:NF]

MIP_SIPM = 0.458
MIP_PL20, MIP_PL25 = 4.334 * 2.0 / 2.5, 4.334   # plastic 0.5-MIP = 1.73 / 2.17 MeV
N_PULSE = 4.284e6
A = 0.5   # SiPM threshold [MIP]
B = 0.5   # plastic threshold [MIP]

tl20 = tl25 = ts = tp20 = tp25 = 0; nev = 0
for fp in FILES:
    obs, aux, n = digest_file(fp, gate=True)      # in-gate (>1 ms)
    nev += n
    S   = obs[:, :, 0] / MIP_SIPM                  # (N,4) per-arm max SiPM [MIP]
    P20 = obs[:, :, 1] / MIP_PL20
    P25 = obs[:, :, 1] / MIP_PL25
    sa  = S >= A
    tl20 += int((sa & (P20 >= B)).sum())           # legs @2.0 cm
    tl25 += int((sa & (P25 >= B)).sum())           # legs @2.5 cm
    ts   += int(sa.sum())                          # SiPM arm-singles
    tp20 += int((P20 >= B).sum())                  # plastic arm-singles @2.0
    tp25 += int((P25 >= B).sum())

f = N_PULSE / nev
print(f"files={len(FILES)}  events={nev:,}  (in-gate, SiPM>={A}, plastic>={B} MIP)")
print(f"LEGS (SiPM^plastic per arm) @2.0cm : {tl20*f:,.0f} /pulse   ({tl20*f/4:,.0f}/arm)")
print(f"LEGS (SiPM^plastic per arm) @2.5cm : {tl25*f:,.0f} /pulse   ({tl25*f/4:,.0f}/arm)")
print(f"  SiPM  arm-singles                : {ts*f:,.0f} /pulse")
print(f"  plastic arm-singles @2.0cm       : {tp20*f:,.0f} /pulse")
print(f"  plastic arm-singles @2.5cm       : {tp25*f:,.0f} /pulse")
