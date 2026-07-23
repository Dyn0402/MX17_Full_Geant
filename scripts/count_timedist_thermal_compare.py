#!/usr/bin/env python3
"""count_timedist_thermal_compare.py — batch-digest count_timedist_thermal.py
over the baseline + no-Al + gamma-cut-scan campaigns, so each writes its own
npz under analysis/al_pair_crosscheck/. Run on lxplus after the campaigns
land (source scripts/setup_lxplus.sh first).

Usage:
    python3 scripts/count_timedist_thermal_compare.py [NF]
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BASE = "/eos/experiment/ntof/data/x17/full_sim"
OUT  = ROOT / "analysis/al_pair_crosscheck"
OUT.mkdir(parents=True, exist_ok=True)

NF = sys.argv[1] if len(sys.argv) > 1 else "8"

RUNS = [
    ("baseline_100um",  f"{BASE}/neutrons_thermal_trig_2cm"),
    ("noAl",            f"{BASE}/neutrons_thermal_trig_2cm_noAl"),
    ("gcut50um",        f"{BASE}/neutrons_thermal_trig_2cm_gcut50um"),
    ("gcut20um",        f"{BASE}/neutrons_thermal_trig_2cm_gcut20um"),
    ("gcut10um",        f"{BASE}/neutrons_thermal_trig_2cm_gcut10um"),
    ("gcut5um",         f"{BASE}/neutrons_thermal_trig_2cm_gcut5um"),
]

for name, d in RUNS:
    out = OUT / f"timedist_{name}.npz"
    print(f"=== {name}  ({d}) -> {out} ===", flush=True)
    subprocess.run(
        [sys.executable, str(HERE / "count_timedist_thermal.py"), d, NF, str(out)],
        check=True,
    )
print("done.")
