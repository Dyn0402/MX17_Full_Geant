#!/usr/bin/env python3
"""
submit_ill.py — HTCondor submission for the ILL campaign (HANDOFF_SIM.md, in
x17_facility_search/ill).  One call submits one (run, configuration) pair.

Configurations G1–G6 are the §4 scan (pressure × wall × radius, each cell
stopping 99.5 % of the H113 beam, entrance window at y_w = −median depth).
Common to all: Be 0.5 mm window, Al 8 mm end cap, ⁶LiF 12:5 scraper, the
realistic PF1B beam (--beam ill data/beam/h113_spectrum.csv).  Anything after
`--` is appended to the simulation arguments (e.g. a C2 window variant, or
`--cell-yw 0` for the S1p placement scan).

    python3 scripts/submit_ill.py --run V0 --config G1 --njobs 10 --nevents 100000
    python3 scripts/submit_ill.py --run C1 --config G3 --njobs 10 --nevents 10000000 \\
        --bias-ncapture 1e5
    python3 scripts/submit_ill.py --run C2 --config G1 --tag win_Al0.1 -- --window Al:0.1

Outputs: /eos/experiment/ntof/data/x17/ill/<run>/<config>[_<tag>]/<prefix>_jobNNN_t0.root
Condor files + logs: /afs/cern.ch/user/d/dneff/condor/ill/<run>/<config>[_<tag>]/
"""

import argparse
import os
import random
import stat
import sys
import textwrap
import zlib
from pathlib import Path

EOS_BASE = "/eos/experiment/ntof/data/x17/ill"
JOB_BASE = "/afs/cern.ch/user/d/dneff/condor/ill"
REPO = Path(__file__).resolve().parent.parent
SPECTRUM = REPO / "data" / "beam" / "h113_spectrum.csv"

COMMON = ["--window", "Be:0.5", "--end-cap", "Al:8", "--scraper", "12:5"]

# id: (p_bar, L_mm, R_mm, y_w_mm, skin, rods)   — HANDOFF_SIM.md §4
CONFIGS = {
    "G1": (1, 300, 40,  -21.6, "Mylar:0.012", 6),
    "G2": (1, 300, 100, -21.6, "Mylar:0.012", 15),
    "G3": (2, 150, 40,  -10.8, "Kapton:0.11", 0),
    "G4": (2, 150, 100, -10.8, "Kapton:0.29", 0),
    "G5": (3, 100, 40,  -7.2,  "Kapton:0.23", 0),
    "G6": (3, 100, 100, -7.2,  "Kapton:0.57", 0),
}


def config_args(cfg):
    if cfg == "capsule":
        return ["--target", "capsule"]
    p, L, R, yw, skin, rods = CONFIGS[cfg]
    a = ["--target", "cell", "--cell-pressure", f"{p:g}", "--cell-length", f"{L:g}",
         "--cell-radius", f"{R:g}", "--cell-yw", f"{yw:g}", "--skin", skin] + COMMON
    if rods:
        a += ["--rods", str(rods)]
    return a


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--run", required=True, help="campaign step: V0, C1, C2, S1, ...")
    ap.add_argument("--config", required=True, help="G1..G6 or capsule")
    ap.add_argument("--tag", default="", help="suffix for variants (dir name)")
    ap.add_argument("--mode", choices=["neutron", "pairs"], default="neutron")
    ap.add_argument("--njobs", type=int, default=10)
    ap.add_argument("--nevents", type=int, default=10_000_000)
    ap.add_argument("--bias-ncapture", type=float, default=1.0)
    ap.add_argument("--flavour", default="workday")
    ap.add_argument("--seed", type=int, default=None, help="master seed (default: crc32 of run/config[_tag])")
    ap.add_argument("--memory", type=int, default=2048)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("extra", nargs=argparse.REMAINDER, help="-- extra simulation args")
    a = ap.parse_args()

    extra = a.extra[1:] if a.extra[:1] == ["--"] else a.extra
    name = a.config + (f"_{a.tag}" if a.tag else "")
    outdir = Path(EOS_BASE) / a.run / name
    jobdir = Path(JOB_BASE) / a.run / name
    exe = REPO / "build" / "mx17_full_sim"
    if not exe.is_file():
        sys.exit(f"ERROR: build first ({exe})")

    sim = config_args(a.config)
    if a.mode == "neutron":
        sim = ["--beam", "ill", str(SPECTRUM)] + sim
        if a.bias_ncapture > 1:
            sim += ["--bias-ncapture", f"{a.bias_ncapture:g}"]
    sim += extra
    prefix = "neutrons" if a.mode == "neutron" else "pairs"

    seed0 = a.seed if a.seed is not None else zlib.crc32(f"{a.run}/{name}".encode())
    rng = random.Random(seed0)
    jobs = [(str(outdir / f"{prefix}_job{i:03d}"), rng.randint(1, 2**31 - 1), f"job{i:03d}")
            for i in range(a.njobs)]

    print(f"{a.run}/{name}: {a.njobs} × {a.nevents:,} {a.mode}  → {outdir}")
    print("  args: " + " ".join(sim))
    if a.dry_run:
        return

    outdir.mkdir(parents=True, exist_ok=True)
    (jobdir / "logs").mkdir(parents=True, exist_ok=True)
    wrapper = jobdir / "run.sh"
    simargs = " ".join(f'"{x}"' for x in sim)
    wrapper.write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        set -eo pipefail
        set +u
        source "{REPO}/scripts/setup_lxplus.sh" > /dev/null
        set -u
        OUT="$1"; SEED="$2"
        echo "node $(hostname)  out $OUT  seed $SEED  $(date)"
        "{exe}" -t 1 -n {a.nevents} -o "$OUT" -s "$SEED" {simargs} > "$OUT.log" 2>&1
        tail -4 "$OUT.log"
        echo "done $(date)"
    """))
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC)
    sub = jobdir / "jobs.sub"
    lines = [
        f"executable = {wrapper}",
        f"output     = {jobdir}/logs/$(tag).out",
        f"error      = {jobdir}/logs/$(tag).err",
        f"log        = {jobdir}/logs/condor.log",
        f'+JobFlavour = "{a.flavour}"',
        "request_cpus = 1",
        f"request_memory = {a.memory}",
        'requirements = (OpSysAndVer =?= "AlmaLinux9")',
        "should_transfer_files = NO",
        "arguments = $(outfile) $(seed)",
        "queue outfile,seed,tag from (",
    ] + [f"  {o}, {s}, {t}" for o, s, t in jobs] + [")"]
    sub.write_text("\n".join(lines) + "\n")
    (jobdir / "README").write_text(f"{a.run}/{name}\n{a.njobs} x {a.nevents} {a.mode}\nargs: {' '.join(sim)}\n")
    (outdir / "RUN_INFO.txt").write_text(f"{a.run}/{name}\n{a.njobs} x {a.nevents} {a.mode}\nmaster seed {seed0}\nargs: {' '.join(sim)}\n")
    if os.system(f"condor_submit {sub}") != 0:
        sys.exit("condor_submit failed")


if __name__ == "__main__":
    main()
