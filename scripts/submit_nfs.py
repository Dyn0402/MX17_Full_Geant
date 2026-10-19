#!/usr/bin/env python3
"""
submit_nfs.py — HTCondor submission for the GANIL/NFS fast-neutron runs
(x17_facility_search/ganil_nfs/FEASIBILITY_SIM.md).  One call submits one (run, sample).

Every sample uses the NFS baseline cell (PLAN_GANIL_NFS.md §4.4 build C: 30 bar ³He,
L = 400 mm centred on the origin, R = 50 mm, Be 0.1 mm + CFRP 0.3 mm skin, Be 3 mm
window in a 45 mm aperture, Al 4 mm end cap; pair vertices within the r = 40 mm
beam spot) unless overridden after `--`.  Samples:

    X17_En<E>        X17 (m 16.8) at W = S_n + 0.75 E_n, ⁴He* boosted by β_cm(E_n)   (N1)
    <E1|M1>_En<E>    Born IPC of that multipole at the same W, opening angle ≥ 40°   (N2)
    beam             NFS d(40)+Be neutrons, 0.1–45 MeV, r = 40 mm disc, vacuum tube  (N3)
    cosmic           cosmic μ                                                          (N4)

    python3 scripts/submit_nfs.py --run N1 --sample X17_En2 --njobs 5 --nevents 100000 --exe bin/...
    python3 scripts/submit_nfs.py --run N2 --sample E1_En2 --tag bigP -- --big-plastic 75 75 5 --no-ls --sipm-readout 20 0

Outputs: /eos/experiment/ntof/data/x17/nfs/<run>/<sample>[_<tag>]/<prefix>_jobNNN_t0.root
Condor files + logs: /afs/cern.ch/user/d/dneff/condor/nfs/<run>/<sample>[_<tag>]/
Jobs write straight to EOS (should_transfer_files = NO), so nothing lands on AFS.
"""

import argparse
import os
import random
import re
import stat
import sys
import textwrap
import zlib
from pathlib import Path

EOS_BASE = "/eos/experiment/ntof/data/x17/nfs"
ANALYSIS = "/afs/cern.ch/work/d/dneff/git/x17_nfs/analysis"     # nfs_reduce.py, lnl_time.py, seg_reduce.py
LCG = "/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh"
JOB_BASE = "/afs/cern.ch/user/d/dneff/condor/nfs"
REPO = Path(__file__).resolve().parent.parent

CELL = ["--target", "cell", "--cell-pressure", "30", "--cell-length", "400", "--cell-radius", "50",
        "--cell-yw", "-200", "--skin", "Be:0.1+CFRP:0.3", "--window", "Be:3", "--aperture", "45",
        "--end-cap", "Al:4", "--scraper", "none"]
BASELINE = CELL + ["--pair-r", "40"]
BEAM = ["--beam", str(REPO / "data/beam/nfs_dbe40_gun.csv"), "--beam-radius", "40", "--gun-dist", "1000",
        "--flight-tube", "Vac:55", "--tube-wall", "Al:1", "--tube-window", "Al:0.5"]
BEAM.insert(1, "nfs")
S_N, M_N, M_HE3 = 20.577, 939.56542, 2808.39161     # MeV


def w_beta(e_n):
    """⁴He* excitation (= transition energy to the g.s.) and its velocity along the beam."""
    m_he4 = 3727.379
    s = (M_N + M_HE3) ** 2 + 2 * M_HE3 * e_n
    w = s ** 0.5 - m_he4
    p = (e_n * (e_n + 2 * M_N)) ** 0.5
    return w, p / (e_n + M_N + M_HE3)


def sample_args(sample):
    if m := re.fullmatch(r"X17_En([\d.]+)(?:_m([\d.]+))?", sample):
        w, b = w_beta(float(m.group(1)))
        return ["--energy", f"{w:.4f}", "--mass", m.group(2) or "16.8", "--ipc", "0",
                "--cm-beta", f"{b:.5f}"], "pairs"
    if m := re.fullmatch(r"(M1|E1)_En([\d.]+)", sample):
        w, b = w_beta(float(m.group(2)))
        return ["--energy", f"{w:.4f}", "--ipc", "1", "--ipc-multipole", m.group(1),
                "--ipc-min-open", "40", "--cm-beta", f"{b:.5f}"], "pairs"
    if sample == "beam":
        return BEAM, "neutrons"
    if sample == "cosmic":
        return ["--cosmic"], "cosmics"
    sys.exit(f"unknown sample '{sample}'")


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--run", required=True, help="N1, N2, ...")
    ap.add_argument("--sample", required=True, help="X17_En2, E1_En2, beam, cosmic, ...")
    ap.add_argument("--tag", default="", help="suffix for variants (dir name)")
    ap.add_argument("--njobs", type=int, default=10)
    ap.add_argument("--nevents", type=int, default=100_000)
    ap.add_argument("--flavour", default="workday")
    ap.add_argument("--seed", type=int, default=None, help="master seed (default: crc32 of run/sample[_tag])")
    ap.add_argument("--memory", type=int, default=2048)
    ap.add_argument("--exe", default=None,
                    help="binary to run (default build/mx17_full_sim); use a frozen copy in "
                         "bin/ so a rebuild cannot touch running jobs")
    ap.add_argument("--reduce-in-job", action="store_true",
                    help="simulate into the condor scratch dir, run nfs_reduce.py + lnl_time.py there "
                         "(parts to <outdir>/parts_nfs, parts_time) and delete the ROOT file: "
                         "for the large --force-beam samples (~1.6 kB/event)")
    ap.add_argument("--keep-root", type=int, default=5,
                    help="with --reduce-in-job: the first N jobs also copy their ROOT file to EOS")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("extra", nargs=argparse.REMAINDER, help="-- extra simulation args")
    a = ap.parse_args()

    extra = a.extra[1:] if a.extra[:1] == ["--"] else a.extra
    name = a.sample + (f"_{a.tag}" if a.tag else "")
    outdir = Path(EOS_BASE) / a.run / name
    jobdir = Path(JOB_BASE) / a.run / name
    exe = Path(a.exe).resolve() if a.exe else REPO / "build" / "mx17_full_sim"
    if not exe.is_file():
        sys.exit(f"ERROR: build first ({exe})")

    gen, prefix = sample_args(a.sample)
    sim = (CELL if prefix == "neutrons" else BASELINE) + gen + extra

    seed0 = a.seed if a.seed is not None else zlib.crc32(f"{a.run}/{name}".encode())
    rng = random.Random(seed0)
    jobs = [(str(outdir / f"{prefix}_job{i:03d}"), rng.randint(1, 2**31 - 1), f"job{i:03d}")
            for i in range(a.njobs)]

    print(f"{a.run}/{name}: {a.njobs} × {a.nevents:,} {prefix}  → {outdir}")
    print("  args: " + " ".join(sim))
    if a.dry_run:
        return

    outdir.mkdir(parents=True, exist_ok=True)
    (jobdir / "logs").mkdir(parents=True, exist_ok=True)
    wrapper = jobdir / "run.sh"
    simargs = " ".join(f'"{x}"' for x in sim)
    if a.reduce_in_job:
        run = textwrap.dedent(f"""\
            D=$(dirname "$OUT"); S=$(basename "$OUT")
            W="${{_CONDOR_SCRATCH_DIR:-/tmp/$USER}}/$S"
            "{exe}" -t 1 -n {a.nevents} -o "$W" -s "$SEED" {simargs} > "$OUT.log" 2>&1
            tail -4 "$OUT.log"
            mkdir -p "$D/parts_nfs" "$D/parts_time"
            ( set +u; source {LCG} > /dev/null
              python3 {ANALYSIS}/nfs_reduce.py "${{W}}_t0.root" -o "$D/parts_nfs/${{S}}_t0.npz"
              python3 {ANALYSIS}/lnl_time.py "${{W}}_t0.root" -o "$D/parts_time/${{S}}_t0.npz" )
            if [ "$KEEP" = 1 ]; then cp "${{W}}_t0.root" "$D/"; fi
            rm -f "$W"*
        """)
    else:
        run = textwrap.dedent(f"""\
            "{exe}" -t 1 -n {a.nevents} -o "$OUT" -s "$SEED" {simargs} > "$OUT.log" 2>&1
            tail -4 "$OUT.log"
        """)
    wrapper.write_text(textwrap.dedent("""\
        #!/usr/bin/env bash
        set -eo pipefail
        set +u
        source "{repo}/scripts/setup_lxplus.sh" > /dev/null
        set -u
        OUT="$1"; SEED="$2"; KEEP="${{3:-0}}"
        echo "node $(hostname)  out $OUT  seed $SEED  $(date)"
        """).format(repo=REPO) + run + 'echo "done $(date)"\n')
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
        "arguments = $(outfile) $(seed) $(keep)",
        "queue outfile,seed,tag,keep from (",
    ] + [f"  {o}, {s}, {t}, {int(a.reduce_in_job and i < a.keep_root)}"
         for i, (o, s, t) in enumerate(jobs)] + [")"]
    sub.write_text("\n".join(lines) + "\n")
    info = (f"{a.run}/{name}\n{a.njobs} x {a.nevents} {prefix}\nmaster seed {seed0}\nexe {exe}\nargs: {' '.join(sim)}\n"
            + ("reduced in the job (parts_nfs/, parts_time/); ROOT kept for the first "
               f"{a.keep_root} jobs\n" if a.reduce_in_job else ""))
    (jobdir / "README").write_text(info)
    (outdir / "RUN_INFO.txt").write_text(info)
    if os.system(f"condor_submit {sub}") != 0:
        sys.exit("condor_submit failed")


if __name__ == "__main__":
    main()
