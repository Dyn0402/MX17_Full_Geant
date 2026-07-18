#!/usr/bin/env python3
"""
submit_pairs.py
Submit X17 + IPC pair simulation jobs to HTCondor on lxplus.

Runs mx17_full_sim in pair mode (default: 50% X17, 50% IPC).
Default: 100 jobs × 100k events = 10M events total.

Workflow:
    1. Build:   bash scripts/build.sh
    2. Submit:  python3 scripts/submit_pairs.py [--dry-run]
    3. Monitor: condor_q
    4. Merge:   hadd x17_ipc_merged.root <outdir>/x17_ipc_pairs_job*.root

Output per job: <outdir>/x17_ipc_pairs_job<NNN>.root
Each ROOT file contains HitTree (per-step hits) and EventTree (per-event truth).
EventTree has event_type branch: 0=X17, 1=IPC, for splitting the pools in Python.
"""

import argparse
import os
import random
import stat
import sys
import textwrap
from pathlib import Path

# ── Defaults ──────────────────────────────────────────────────────────────────
OUTDIR_DEFAULT  = "/eos/experiment/ntof/data/x17/full_sim/pairs"
JOBDIR_DEFAULT  = "/afs/cern.ch/user/d/dneff/condor/mx17_pairs"
NJOBS_DEFAULT   = 100
NEVENTS_DEFAULT = 100_000
IPC_DEFAULT     = 0.5      # 50/50 X17 vs IPC
GAS_DEFAULT     = "ArIso"
JOB_FLAVOUR     = "workday"  # ~8 h; increase to "tomorrow" if jobs time out


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Submit X17+IPC pair sim jobs to HTCondor on lxplus",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--dry-run",  action="store_true",
                   help="Print job list without submitting")
    p.add_argument("--njobs",    type=int,   default=NJOBS_DEFAULT,
                   help="Number of Condor jobs")
    p.add_argument("--nevents",  type=int,   default=NEVENTS_DEFAULT,
                   help="Events per job")
    p.add_argument("--ipc",      type=float, default=IPC_DEFAULT,
                   help="IPC fraction: 0 = all X17, 1 = all IPC")
    p.add_argument("--vertex-lib", default=None, metavar="GAS_CSV",
                   help="He3Gas capture-position library (make_capture_library.py "
                        "--gas-lib): sample pair vertices from the thermal "
                        "self-shielding profile instead of uniformly in the gas")
    p.add_argument("--gas",      default=GAS_DEFAULT,
                   help="MM gas mixture")
    p.add_argument("--outdir",   default=OUTDIR_DEFAULT,
                   help="Output directory for ROOT files (EOS path)")
    p.add_argument("--jobdir",   default=JOBDIR_DEFAULT,
                   help="Directory for Condor submit files and logs (AFS path)")
    p.add_argument("--exe",      default=None,
                   help="Path to mx17_full_sim (auto-detected from build/ if omitted)")
    p.add_argument("--flavour",  default=JOB_FLAVOUR,
                   help="HTCondor +JobFlavour")
    p.add_argument("--seed",     type=int,   default=42,
                   help="Master RNG seed for deriving per-job seeds")
    return p.parse_args()


# ── Executable search ─────────────────────────────────────────────────────────
def find_exe():
    script_dir = Path(__file__).parent
    candidates = [
        script_dir.parent / "build"   / "mx17_full_sim",
        script_dir.parent / "install" / "bin" / "mx17_full_sim",
        Path(os.environ.get("MX17_EXE", "")),
    ]
    for c in candidates:
        if c.is_file():
            return str(c.resolve())
    return None


# ── Bash wrapper (one file, shared by all jobs) ───────────────────────────────
def write_wrapper(job_dir: Path, exe: str, setup_script: str,
                  vertex_lib=None) -> Path:
    extra = f'\\\n            --pair-vertex-lib "{vertex_lib}" ' if vertex_lib else ""
    wrapper = job_dir / "run_pairs_job.sh"
    wrapper.write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        set -eo pipefail
        # Disable -u around setup sourcing: LCG/GCC scripts may reference
        # unset variables internally before they define them.
        set +u
        source "{setup_script}"
        set -u

        OUTFILE="$1"
        NEVENTS="$2"
        SEED="$3"
        IPC="$4"
        GAS="$5"

        echo "=== mx17_full_sim pair job ==="
        echo "  Node    : $(hostname)"
        echo "  Output  : $OUTFILE"
        echo "  Events  : $NEVENTS"
        echo "  IPC     : $IPC"
        echo "  Gas     : $GAS"
        echo "  Seed    : $SEED"
        echo "=============================="

        "{exe}" \\
            -n "$NEVENTS" \\
            -o "$OUTFILE" \\
            -s "$SEED"    \\
            --ipc "$IPC"  \\
            -g "$GAS" {extra}

        echo "Job done: $(date)"
    """))
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return wrapper


# ── Condor submit file ────────────────────────────────────────────────────────
def write_submit(job_dir: Path, wrapper: Path, jobs: list, flavour: str) -> Path:
    log_dir = job_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    sub = job_dir / "pairs.sub"

    lines = [
        f"executable             = {wrapper}",
        f"output                 = {log_dir}/$(tag).out",
        f"error                  = {log_dir}/$(tag).err",
        f"log                    = {log_dir}/condor.log",
        f'+JobFlavour            = "{flavour}"',
        "request_cpus           = 1",
        "request_memory         = 2048",
        "request_disk           = 4096",
        "should_transfer_files  = YES",
        "when_to_transfer_output = ON_EXIT",
        "",
        'requirements           = (OpSysAndVer =?= "AlmaLinux9")',
        "",
        "# Positional args passed to wrapper: outfile nevents seed ipc_fraction gas",
        "arguments              = $(outfile) $(nevents) $(seed) $(ipc) $(gas)",
        "",
        "queue outfile,nevents,seed,ipc,gas,tag from (",
    ]
    for (outfile, nevents, seed, ipc, gas, tag) in jobs:
        lines.append(f"  {outfile}, {nevents}, {seed}, {ipc:.3f}, {gas}, {tag}")
    lines.append(")")

    sub.write_text("\n".join(lines) + "\n")
    return sub


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    exe = args.exe or find_exe()
    if not exe or not Path(exe).is_file():
        print("ERROR: mx17_full_sim not found.")
        print("  Build first:  bash scripts/build.sh")
        print("  Or specify:   --exe /path/to/mx17_full_sim")
        sys.exit(1)
    exe = str(Path(exe).resolve())

    outdir  = Path(args.outdir)
    job_dir = Path(args.jobdir)
    outdir.mkdir(parents=True, exist_ok=True)
    job_dir.mkdir(parents=True, exist_ok=True)

    setup_script = str((Path(__file__).parent / "setup_lxplus.sh").resolve())

    rng = random.Random(args.seed)
    jobs = []
    for i in range(args.njobs):
        tag     = f"x17_ipc_pairs_job{i:03d}"
        outfile = str(outdir / tag)
        seed    = rng.randint(1, 2**31 - 1)
        jobs.append((outfile, args.nevents, seed, args.ipc, args.gas, tag))

    total_events = args.njobs * args.nevents
    x17_events   = int(total_events * (1 - args.ipc))
    ipc_events   = int(total_events * args.ipc)

    print(f"Jobs           : {args.njobs}")
    print(f"Events/job     : {args.nevents:,}")
    print(f"Total events   : {total_events:,}  "
          f"(~{x17_events:,} X17  +  ~{ipc_events:,} IPC)")
    print(f"IPC fraction   : {args.ipc}")
    print(f"Gas            : {args.gas}")
    print(f"Job flavour    : {args.flavour}")
    print(f"Output dir     : {outdir}")
    print(f"Condor dir     : {job_dir}")
    print(f"Executable     : {exe}")

    if args.dry_run:
        print("\n--- DRY RUN: first 5 jobs ---")
        for (outfile, nevents, seed, ipc, gas, tag) in jobs[:5]:
            print(f"  {tag}  seed={seed}")
        if args.njobs > 5:
            print(f"  ... and {args.njobs - 5} more")
        print(f"\nMerge command (after completion):")
        print(f"  hadd {outdir}/x17_ipc_merged.root {outdir}/x17_ipc_pairs_job*.root")
        return

    vertex_lib = str(Path(args.vertex_lib).resolve()) if args.vertex_lib else None
    if vertex_lib and not Path(vertex_lib).is_file():
        print(f"ERROR: vertex library not found: {vertex_lib}")
        sys.exit(1)
    wrapper  = write_wrapper(job_dir, exe, setup_script, vertex_lib)
    sub_file = write_submit(job_dir, wrapper, jobs, args.flavour)

    print(f"\nSubmit file : {sub_file}")
    print("Submitting to HTCondor...")
    ret = os.system(f"condor_submit {sub_file}")
    if ret != 0:
        print("ERROR: condor_submit failed")
        sys.exit(1)

    print(f"\nDone. Monitor with:  condor_q")
    print(f"\nAfter all jobs complete, merge outputs:")
    print(f"  hadd {outdir}/x17_ipc_merged.root {outdir}/x17_ipc_pairs_job*.root")


if __name__ == "__main__":
    main()
