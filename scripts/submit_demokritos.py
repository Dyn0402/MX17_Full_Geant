#!/usr/bin/env python3
"""
submit_demokritos.py — HTCondor submission for the Demokritos ³H(p,e⁺e⁻)⁴He runs
(x17_facility_search/demokritos/GEANT_PREP.md).  One call submits one (run, sample).

The --target li stack with a Ti-tritide film: TiT₁.₅ 2300 µg/cm² (2.1 mg/cm² Ti,
Demokritos' own targets) on a Cu backing (default 50 µm; --backing Cu:1000 is the
existing target), Al 1 mm holder, CFRP 0.4 mm chamber of bore 25 mm, no beam dump
(the protons stop in the tritide / backing).  Override after `--`.  Samples:

    X17_m<m>         X17 → e⁺e⁻ at W = 20.49 MeV (E_p = 900 keV), mass m   (D1)
    <E1|M1|E0>_<W>   Born IPC of that multipole at transition W            (D2)
    gam_<name>       --gamma-lines; name keys GAMMA_LINES                  (D3)
    cosmic           cosmic μ                                              (reuse LNL L4)

    python3 scripts/submit_demokritos.py --run D1 --sample X17_m16.94 --tag cu50 --njobs 10 \\
        --nevents 100000 --exe bin/mx17_full_sim_demok_<sha>
    python3 scripts/submit_demokritos.py --run D1 --sample X17_m16.94 --tag cu1000_y120 -- \\
        --backing Cu:1000 --target-y 120

Outputs: /eos/experiment/ntof/data/x17/demokritos/<run>/<sample>[_<tag>]/<prefix>_jobNNN_t0.root
Condor files + logs: /afs/cern.ch/user/d/dneff/condor/demokritos/<run>/<sample>[_<tag>]/
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

EOS_BASE = "/eos/experiment/ntof/data/x17/demokritos"
JOB_BASE = "/afs/cern.ch/user/d/dneff/condor/demokritos"
REPO = Path(__file__).resolve().parent.parent

BASELINE = ["--target", "li", "--film", "TiT:2300", "--backing", "Cu:50", "--holder", "Al:1",
            "--chamber", "CFRP:0.4:25", "--dump", "none", "--spot-sigma", "2"]
W_RES = 20.49            # MeV, E* at E_p = 900 keV (E* = 19.814 + 0.75 E_p)

# --gamma-lines sets (E_MeV:weight); weights are only the sampling mix, every
# event carries its E_γ so the analysis reweights to the real yields.
GAMMA_LINES = {
    # the 3H(p,g)4He ground-state line over the thick-target E* range (EPC source)
    "He4": "20.49:1",
    "He4w": "20.2:1,20.35:1,20.49:1",
}


def sample_args(sample):
    if m := re.fullmatch(r"X17_m([\d.]+)", sample):
        return ["--energy", f"{W_RES:g}", "--mass", m.group(1), "--ipc", "0"], "pairs"
    if m := re.fullmatch(r"(M1|E1|E0)_([\d.]+)", sample):
        return ["--energy", m.group(2), "--ipc", "1", "--ipc-multipole", m.group(1)], "pairs"
    if m := re.fullmatch(r"gam_(\w+)", sample):
        return ["--gamma-lines", GAMMA_LINES[m.group(1)]], "gammas"
    if sample == "cosmic":
        return ["--cosmic"], "cosmics"
    sys.exit(f"unknown sample '{sample}'")


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--run", required=True, help="D1, D2, ... (demokritos/GEANT_PREP.md)")
    ap.add_argument("--sample", required=True, help="X17_m16.94, E1_20.49, gam_He4, ...")
    ap.add_argument("--tag", default="", help="suffix for variants (dir name)")
    ap.add_argument("--njobs", type=int, default=10)
    ap.add_argument("--nevents", type=int, default=100_000)
    ap.add_argument("--flavour", default="workday")
    ap.add_argument("--seed", type=int, default=None, help="master seed (default: crc32 of run/sample[_tag])")
    ap.add_argument("--memory", type=int, default=2048)
    ap.add_argument("--exe", default=None,
                    help="binary to run (default build/mx17_full_sim); use a frozen copy in "
                         "bin/ so a rebuild cannot touch running jobs")
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
    sim = BASELINE + gen + extra

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
    info = f"{a.run}/{name}\n{a.njobs} x {a.nevents} {prefix}\nmaster seed {seed0}\nexe {exe}\nargs: {' '.join(sim)}\n"
    (jobdir / "README").write_text(info)
    (outdir / "RUN_INFO.txt").write_text(info)
    if os.system(f"condor_submit {sub}") != 0:
        sys.exit("condor_submit failed")


if __name__ == "__main__":
    main()
