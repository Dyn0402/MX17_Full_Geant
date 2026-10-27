#!/usr/bin/env python3
"""
submit_room.py — HTCondor submission for the Birmingham room runs (--room; x17_facility_search
birmingham/TARGET_GEANT.md §3): the source → shield-wall stages and the experimental-room runs.

    python3 scripts/submit_room.py --run R1 --name s1_n --mode ps --njobs 200 --nevents 200000 --exe bin/... \
        -- --room --src data/beam/bham_src_n_28.csv --ps-plane -1450
    python3 scripts/submit_room.py --run R1 --name s2_n --mode ps ... -- --room --ps-gun <merged s1 csv> --ps-plane ...
    python3 scripts/submit_room.py --run R2 --name leak_n_bigP --mode det ... -- --room --ps-gun ... <detector args>

--mode ps   each job copies its <out>_ps_t0.csv to EOS (bham_ps_merge.py then normalises the stage)
--mode det  each job reduces in the scratch dir (nfs_beam_reduce.py → parts_beam, bham_chan_reduce.py →
            parts_chan) and keeps the ROOT file for the first --keep-root jobs
Every finished job writes <out>.nev (events simulated), so the merge knows the denominator.
Outputs /eos/experiment/ntof/data/x17/bham/<run>/<name>/ ; condor files /afs/cern.ch/user/d/dneff/condor/bham/<run>/<name>/
"""

import argparse
import os
import random
import stat
import sys
import textwrap
import zlib
from pathlib import Path

EOS_BASE = "/eos/experiment/ntof/data/x17/bham"
ANALYSIS = "/afs/cern.ch/work/d/dneff/git/x17_nfs/analysis"
LCG = "/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh"
JOB_BASE = "/afs/cern.ch/user/d/dneff/condor/bham"
REPO = Path(__file__).resolve().parent.parent


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--mode", choices=("ps", "det"), required=True)
    ap.add_argument("--njobs", type=int, default=100)
    ap.add_argument("--nevents", type=int, default=100_000)
    ap.add_argument("--flavour", default="workday")
    ap.add_argument("--memory", type=int, default=2048)
    ap.add_argument("--exe", required=True)
    ap.add_argument("--keep-root", type=int, default=2)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("extra", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    sim = a.extra[1:] if a.extra[:1] == ["--"] else a.extra
    outdir, jobdir = Path(EOS_BASE) / a.run / a.name, Path(JOB_BASE) / a.run / a.name
    exe = Path(a.exe).resolve()
    if not exe.is_file():
        sys.exit(f"ERROR: no binary {exe}")
    seed0 = zlib.crc32(f"{a.run}/{a.name}".encode())
    rng = random.Random(seed0)
    jobs = [(str(outdir / f"room_job{i:03d}"), rng.randint(1, 2**31 - 1), f"job{i:03d}") for i in range(a.njobs)]
    print(f"{a.run}/{a.name}: {a.njobs} × {a.nevents:,} ({a.mode}) → {outdir}\n  args: {' '.join(sim)}")
    if a.dry_run:
        return
    outdir.mkdir(parents=True, exist_ok=True)
    (jobdir / "logs").mkdir(parents=True, exist_ok=True)
    simargs = " ".join(f'"{x}"' for x in sim)
    if a.mode == "ps":
        post = 'cp "${W}_ps_t0.csv" "$OUT.ps.csv" 2>/dev/null || true\n'
    else:
        post = textwrap.dedent(f"""\
            mkdir -p "$D/parts_beam" "$D/parts_chan"
            ( set +u; source {LCG} > /dev/null
              python3 {ANALYSIS}/nfs_beam_reduce.py "${{W}}_t0.root" -o "$D/parts_beam/${{S}}_t0.npz"
              python3 {ANALYSIS}/bham_chan_reduce.py "${{W}}_t0.root" -o "$D/parts_chan/${{S}}_t0.npz" )
            if [ "$KEEP" = 1 ]; then cp "${{W}}_t0.root" "$D/"; fi
        """)
    wrapper = jobdir / "run.sh"
    wrapper.write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        set -eo pipefail
        set +u
        source "{REPO}/scripts/setup_lxplus.sh" > /dev/null
        set -u
        OUT="$1"; SEED="$2"; KEEP="${{3:-0}}"
        echo "node $(hostname)  out $OUT  seed $SEED  $(date)"
        D=$(dirname "$OUT"); S=$(basename "$OUT")
        W="${{_CONDOR_SCRATCH_DIR:-/tmp/$USER}}/$S"
        cd "$(dirname "$W")"
        "{exe}" -t 1 -n {a.nevents} -o "$W" -s "$SEED" {simargs} > "$OUT.log" 2>&1
        tail -2 "$OUT.log"
        """) + post + textwrap.dedent(f"""\
        echo {a.nevents} > "$OUT.nev"
        rm -f "$W"*
        echo "done $(date)"
        """))
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC)
    sub = jobdir / "jobs.sub"
    sub.write_text("\n".join([
        f"executable = {wrapper}",
        f"output     = {jobdir}/logs/$(tag).out",
        f"error      = {jobdir}/logs/$(tag).err",
        f"log        = {jobdir}/logs/condor.log",
        f'+JobFlavour = "{a.flavour}"',
        "request_cpus = 1",
        f"request_memory = {a.memory}",
        'requirements = (OpSysAndVer =?= "AlmaLinux9")',
        "should_transfer_files = NO",
        "arguments = $(outfile) $(seed) $(keep)",
        "queue outfile,seed,tag,keep from (",
    ] + [f"  {o}, {s}, {t}, {int(i < a.keep_root)}" for i, (o, s, t) in enumerate(jobs)] + [")"]) + "\n")
    info = (f"{a.run}/{a.name}\n{a.njobs} x {a.nevents} {a.mode}\nmaster seed {seed0}\nexe {exe}\n"
            f"args: {' '.join(sim)}\n")
    (jobdir / "README").write_text(info)
    (outdir / "RUN_INFO.txt").write_text(info)
    for _ in range(4):
        if os.system(f"condor_submit {sub}") == 0:
            return
    sys.exit("condor_submit failed")


if __name__ == "__main__":
    main()
