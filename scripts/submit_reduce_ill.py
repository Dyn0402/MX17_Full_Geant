#!/usr/bin/env python3
"""
submit_reduce_ill.py — one condor job per ROOT file of an ILL run directory,
running ill_accounting.py reduce (neutron runs) or ill_pairs.py reduce (pair
runs) into <dir>/parts/.  Files whose part already exists are skipped.

    python3 scripts/submit_reduce_ill.py /eos/experiment/ntof/data/x17/ill/C1/G1 --kind accounting
    python3 scripts/submit_reduce_ill.py /eos/experiment/ntof/data/x17/ill/S1/G1_X17 --kind pairs
"""
import argparse
import os
import stat
import sys
import textwrap
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
JOB_BASE = Path("/afs/cern.ch/user/d/dneff/condor/ill/reduce")
LCG = "/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--kind", choices=["accounting", "pairs"], required=True)
    ap.add_argument("--flavour", default="longlunch")
    ap.add_argument("--memory", type=int, default=4000)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    script = REPO / "scripts" / ("ill_accounting.py" if a.kind == "accounting" else "ill_pairs.py")
    ext = ".json" if a.kind == "accounting" else ".npz"
    for d in map(Path, a.dirs):
        files = sorted(d.glob("*_t0.root"))
        parts = d / "parts"
        todo = [f for f in files if not (parts / (f.stem + ext)).exists()]
        name = "_".join(d.parts[-2:])
        print(f"{d}: {len(files)} files, {len(todo)} to reduce")
        if not todo or a.dry_run:
            continue
        parts.mkdir(exist_ok=True)
        jd = JOB_BASE / name
        (jd / "logs").mkdir(parents=True, exist_ok=True)
        w = jd / "run.sh"
        w.write_text(textwrap.dedent(f"""\
            #!/usr/bin/env bash
            set -e
            source {LCG} > /dev/null
            python3 {script} reduce "$1" -o "$2"
        """))
        w.chmod(w.stat().st_mode | stat.S_IEXEC)
        sub = jd / "reduce.sub"
        lines = [f"executable = {w}",
                 f"output = {jd}/logs/$(tag).out", f"error = {jd}/logs/$(tag).err",
                 f"log = {jd}/logs/condor.log", f'+JobFlavour = "{a.flavour}"',
                 f"request_memory = {a.memory}", 'requirements = (OpSysAndVer =?= "AlmaLinux9")',
                 "should_transfer_files = NO", "arguments = $(inp) $(out)",
                 "queue inp,out,tag from ("]
        lines += [f"  {f}, {parts / (f.stem + ext)}, {f.stem}" for f in todo] + [")"]
        sub.write_text("\n".join(lines) + "\n")
        if os.system(f"condor_submit {sub}") != 0:
            sys.exit("condor_submit failed")


if __name__ == "__main__":
    main()
