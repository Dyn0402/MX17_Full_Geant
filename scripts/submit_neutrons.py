#!/usr/bin/env python3
"""
submit_neutrons.py
Submit neutron-beam background jobs (and optional gamma-source jobs) to
HTCondor on lxplus.

Run B (default):    mx17_full_sim --neutron <flux> <profile> --emin --emax
    100 jobs × 10M neutrons = 1e9 neutrons through the realistic EAR2 beam
    (energy from the Ph3 evaluated flux, radius from Lambda2D).
    Output: capture budget + singles + unbiased wall-conversion background.

Run C (--gamma-source <lib.csv>):  biased wall-background gammas
    Requires a capture library first:
        python3 scripts/make_capture_library.py <runB outputs> -o capture_lib.csv

Workflow:
    1. Build:    bash scripts/build.sh
    2. Run B:    python3 scripts/submit_neutrons.py [--dry-run]
    3. Library:  python3 scripts/make_capture_library.py \
                     /eos/.../neutrons/*.root -o capture_lib.csv
    4. Run C:    python3 scripts/submit_neutrons.py --gamma-source capture_lib.csv
"""

import argparse
import os
import random
import stat
import sys
import textwrap
from pathlib import Path

OUTDIR_DEFAULT  = "/eos/experiment/ntof/data/x17/full_sim/neutrons"
JOBDIR_DEFAULT  = "/afs/cern.ch/user/d/dneff/condor/mx17_neutrons"
NJOBS_DEFAULT   = 100
NEVENTS_DEFAULT = 10_000_000
EMIN_DEFAULT    = 1e-3      # eV
EMAX_DEFAULT    = 1000.0    # eV  (X17 ROI: E_n < 1 keV)
GAS_DEFAULT     = "ArIso"
JOB_FLAVOUR     = "workday"

DATA_DIR     = Path(__file__).resolve().parent.parent / "data"
FLUX_FILE    = DATA_DIR / "fluxEAR2-Ph3_in_different_units.root"
PROFILE_FILE = DATA_DIR / "lamda2DvsEn_EAR2.root"


def parse_args():
    p = argparse.ArgumentParser(
        description="Submit neutron-beam / gamma-source background jobs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--dry-run",  action="store_true")
    p.add_argument("--njobs",    type=int,   default=NJOBS_DEFAULT)
    p.add_argument("--nevents",  type=int,   default=NEVENTS_DEFAULT,
                   help="Primaries per job")
    p.add_argument("--emin",     type=float, default=EMIN_DEFAULT,
                   help="Neutron sampling window min [eV]")
    p.add_argument("--emax",     type=float, default=EMAX_DEFAULT,
                   help="Neutron sampling window max [eV]")
    p.add_argument("--bias-ncapture", type=float, default=1.0, metavar="FACTOR",
                   help="Scale ³He(n,γ) cross-section in the gas by FACTOR "
                        "(variance reduction; each capture weighted 1/FACTOR). "
                        "1.0 = analog/off. Recommend 1e5–1e6.")
    p.add_argument("--gas",      default=GAS_DEFAULT)
    p.add_argument("--gamma-source", default=None, metavar="LIB_CSV",
                   help="Submit gamma-source (run C) instead of neutrons, "
                        "using this capture library")
    p.add_argument("--outdir",   default=OUTDIR_DEFAULT)
    p.add_argument("--jobdir",   default=JOBDIR_DEFAULT)
    p.add_argument("--exe",      default=None)
    p.add_argument("--flavour",  default=JOB_FLAVOUR)
    p.add_argument("--seed",     type=int, default=42)
    return p.parse_args()


def find_exe():
    script_dir = Path(__file__).parent
    for c in [script_dir.parent / "build" / "mx17_full_sim",
              script_dir.parent / "install" / "bin" / "mx17_full_sim",
              Path(os.environ.get("MX17_EXE", ""))]:
        if c.is_file():
            return str(c.resolve())
    return None


def write_wrapper(job_dir: Path, exe: str, setup_script: str,
                  gamma_lib, bias_ncapture: float = 1.0) -> Path:
    """One wrapper for both modes; mode-specific args baked in."""
    if gamma_lib:
        mode_args = f'--gamma-source "{gamma_lib}"'
        name = "run_gamma_job.sh"
    else:
        bias_arg = (f' --bias-ncapture {bias_ncapture:g}'
                    if bias_ncapture > 1.0 else '')
        mode_args = (f'--neutron "{FLUX_FILE}" "{PROFILE_FILE}" '
                     f'--emin "$EMIN" --emax "$EMAX"{bias_arg}')
        name = "run_neutron_job.sh"
    wrapper = job_dir / name
    wrapper.write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        set -eo pipefail
        set +u
        source "{setup_script}"
        set -u

        OUTFILE="$1"
        NEVENTS="$2"
        SEED="$3"
        GAS="$4"
        EMIN="${{5:-0}}"
        EMAX="${{6:-0}}"

        echo "=== mx17_full_sim background job ==="
        echo "  Node    : $(hostname)"
        echo "  Output  : $OUTFILE"
        echo "  Events  : $NEVENTS"
        echo "  Seed    : $SEED"
        echo "===================================="

        "{exe}" \\
            -n "$NEVENTS" \\
            -o "$OUTFILE" \\
            -s "$SEED"    \\
            -g "$GAS"     \\
            {mode_args}

        echo "Job done: $(date)"
    """))
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return wrapper


def write_submit(job_dir: Path, wrapper: Path, jobs: list, flavour: str) -> Path:
    log_dir = job_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    sub = job_dir / "neutrons.sub"
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
        "arguments              = $(outfile) $(nevents) $(seed) $(gas) $(emin) $(emax)",
        "",
        "queue outfile,nevents,seed,gas,emin,emax,tag from (",
    ]
    for (outfile, nevents, seed, gas, emin, emax, tag) in jobs:
        lines.append(f"  {outfile}, {nevents}, {seed}, {gas}, {emin}, {emax}, {tag}")
    lines.append(")")
    sub.write_text("\n".join(lines) + "\n")
    return sub


def main():
    args = parse_args()
    gamma_mode = args.gamma_source is not None

    exe = args.exe or find_exe()
    if not exe or not Path(exe).is_file():
        print("ERROR: mx17_full_sim not found. Build first: bash scripts/build.sh")
        sys.exit(1)
    exe = str(Path(exe).resolve())

    if gamma_mode:
        lib = Path(args.gamma_source).resolve()
        if not lib.is_file():
            print(f"ERROR: capture library not found: {lib}")
            sys.exit(1)
    else:
        for f in (FLUX_FILE, PROFILE_FILE):
            if not f.is_file():
                print(f"ERROR: beam data file not found: {f}")
                print("  (data/ directory must be present in the repo clone)")
                sys.exit(1)

    outdir  = Path(args.outdir)
    job_dir = Path(args.jobdir)
    if not args.dry_run:
        outdir.mkdir(parents=True, exist_ok=True)
        job_dir.mkdir(parents=True, exist_ok=True)
    setup_script = str((Path(__file__).parent / "setup_lxplus.sh").resolve())

    prefix = "wall_gammas" if gamma_mode else "neutron_bg"
    rng = random.Random(args.seed)
    jobs = []
    for i in range(args.njobs):
        tag     = f"{prefix}_job{i:03d}"
        outfile = str(outdir / tag)
        seed    = rng.randint(1, 2**31 - 1)
        jobs.append((outfile, args.nevents, seed, args.gas,
                     args.emin, args.emax, tag))

    print(f"Mode           : {'gamma-source (run C)' if gamma_mode else 'neutron beam (run B)'}")
    print(f"Jobs           : {args.njobs}")
    print(f"Events/job     : {args.nevents:,}")
    print(f"Total primaries: {args.njobs * args.nevents:,}")
    if not gamma_mode:
        print(f"Energy window  : [{args.emin}, {args.emax}] eV")
        if args.bias_ncapture > 1.0:
            print(f"nCapture bias  : ×{args.bias_ncapture:g}  (³He(n,γ) in gas; "
                  f"weight 1/factor)")
        print(f"Flux file      : {FLUX_FILE}")
        print(f"Profile file   : {PROFILE_FILE}")
    else:
        print(f"Capture library: {args.gamma_source}")
    print(f"Output dir     : {outdir}")
    print(f"Executable     : {exe}")

    if args.dry_run:
        print("\n--- DRY RUN: first 3 jobs ---")
        for j in jobs[:3]:
            print(f"  {j[-1]}  seed={j[2]}")
        return

    wrapper  = write_wrapper(job_dir, exe, setup_script, args.gamma_source,
                             args.bias_ncapture)
    sub_file = write_submit(job_dir, wrapper, jobs, args.flavour)
    print(f"\nSubmit file : {sub_file}")
    ret = os.system(f"condor_submit {sub_file}")
    if ret != 0:
        print("ERROR: condor_submit failed")
        sys.exit(1)
    print("\nDone. Monitor with:  condor_q")
    if not gamma_mode:
        print("\nNext steps after completion:")
        print(f"  python3 scripts/make_capture_library.py {outdir}/ -o capture_lib.csv")
        print(f"  python3 scripts/submit_neutrons.py --gamma-source capture_lib.csv")


if __name__ == "__main__":
    main()
