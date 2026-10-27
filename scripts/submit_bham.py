#!/usr/bin/env python3
"""
submit_bham.py — HTCondor submission for the Birmingham HF-ADNeF target study
(x17_facility_search/birmingham/TARGET_GEANT.md).  One call submits one (run, sample, target, detector).

Targets (the ³He in the collimated ⁷Li(p,n) beam; --target is the key):

    cap500   the n_TOF 500 bar capsule as built (default target), Ø20 mm beam, 300 mm of air before it
    c30      NFS build C: 30 bar × 400 mm, R 50, Be 0.1 + CFRP 0.3 skin, Be 3 window, Al 4 end cap, Ø80 beam
    c10      10 bar × 400 mm, R 50, Al 25 µm + CFRP 0.3 skin, Be 3 window, Al 2 end cap, Ø80 beam
    bal5     1 bar mylar balloon, R 50 × 400 mm: PET 12 µm + Al 7 µm laminate on 6 CFRP rods, PET 100 µm
             end windows, CFRP 3 mm end ring; Ø80 beam
    bal10    the same balloon at R 100 (12.6 L(STP), the capsule's ³He); Ø180 beam

Samples:

    X17_En<E>, E1_En<E>, M1_En<E>   pairs at W = E*(E_n), boosted by β_cm (as submit_nfs.py)
    n          neutrons flat 0.01–1.3 MeV (data/beam/bham_flat_n.csv), --force-beam --roulette 0.2 0.05
    g478       478 keV photons (data/beam/bham_478_g.csv), --roulette 0.2 0.05
    g8be       photons flat 14–19 MeV (data/beam/soft_flat14_19_g.csv), --roulette 0.2 0.05

Detectors: asb (as built, 16 SiPM bars read out) and bigP (75 × 75 × 5 cm slab, 20 bars, no LS).

    python3 scripts/submit_bham.py --run B1 --sample X17_En0.5 --target bal10 --det bigP --njobs 10 --nevents 30000 --exe bin/...

Outputs: /eos/experiment/ntof/data/x17/bham/<run>/<sample>_<target>_<det>/ ; reduced in the job
(pairs: parts_nfs + parts_time; beam: parts_beam + parts_chan, bham_chan_reduce.py).
Condor files + logs: /afs/cern.ch/user/d/dneff/condor/bham/<run>/<sample>_<target>_<det>/
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

EOS_BASE = "/eos/experiment/ntof/data/x17/bham"
ANALYSIS = "/afs/cern.ch/work/d/dneff/git/x17_nfs/analysis"
LCG = "/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh"
JOB_BASE = "/afs/cern.ch/user/d/dneff/condor/bham"
REPO = Path(__file__).resolve().parent.parent

TUBE = ["--tube-wall", "Al:1", "--tube-window", "Al:0.5"]
BALLOON = ["--target", "cell", "--cell-pressure", "1", "--cell-length", "400", "--cell-yw", "-200",
           "--skin", "Mylar:0.012+Al:0.007", "--rods", "6", "--window", "Mylar:0.1", "--end-cap", "Mylar:0.1",
           "--ring", "CFRP", "--endcap-ring", "3:6", "--scraper", "none"]
# target -> (geometry args, beam radius mm, gun distance mm, flight tube or None)
TARGETS = {
    "cap500": (["--target", "capsule"], 10, 300, None),
    "c30": (["--target", "cell", "--cell-pressure", "30", "--cell-length", "400", "--cell-radius", "50",
             "--cell-yw", "-200", "--skin", "Be:0.1+CFRP:0.3", "--window", "Be:3", "--aperture", "45",
             "--end-cap", "Al:4", "--scraper", "none"], 40, 1000, "Vac:55"),
    "c10": (["--target", "cell", "--cell-pressure", "10", "--cell-length", "400", "--cell-radius", "50",
             "--cell-yw", "-200", "--skin", "Al:0.025+CFRP:0.3", "--window", "Be:3", "--aperture", "45",
             "--end-cap", "Al:2", "--scraper", "none"], 40, 1000, "Vac:55"),
    "bal5": (BALLOON + ["--cell-radius", "50", "--aperture", "45"], 40, 1000, "Vac:55"),
    "bal10": (BALLOON + ["--cell-radius", "100", "--aperture", "95"], 90, 1000, "Vac:105"),
}
DETS = {"asb": ["--sipm-readout", "16", "1"],
        "bigP": ["--big-plastic", "75", "75", "5", "--no-ls", "--sipm-readout", "20", "0"]}
GUN = {"n": ("data/beam/bham_flat_n.csv", "neutron", True),
       "g478": ("data/beam/bham_478_g.csv", "gamma", False),
       "g8be": ("data/beam/soft_flat14_19_g.csv", "gamma", False)}
S_N, M_N, M_HE3 = 20.577, 939.56542, 2808.39161     # MeV


def w_beta(e_n):
    """⁴He* excitation (= transition energy to the g.s.) and its velocity along the beam."""
    m_he4 = 3727.379
    s = (M_N + M_HE3) ** 2 + 2 * M_HE3 * e_n
    w = s ** 0.5 - m_he4
    p = (e_n * (e_n + 2 * M_N)) ** 0.5
    return w, p / (e_n + M_N + M_HE3)


def target_args(target, beam):
    geo, rb, gun, tube = TARGETS[target]
    out = list(geo)
    if beam:
        out += ["--beam", "nfs", str(REPO / GUN[beam][0]), "--beam-radius", str(rb), "--gun-dist", str(gun)]
        if tube:
            out += ["--flight-tube", tube] + TUBE
        if GUN[beam][1] == "gamma":
            out += ["--beam-particle", "gamma", "--roulette", "0.2", "0.05"]
        else:
            out += ["--force-beam", "--roulette", "0.2", "0.05"]
    elif geo[1] == "cell":
        out += ["--pair-r", str(rb)]
    return out


def sample_args(sample):
    if m := re.fullmatch(r"X17_En([\d.]+)(?:_m([\d.]+))?", sample):
        w, b = w_beta(float(m.group(1)))
        return ["--energy", f"{w:.4f}", "--mass", m.group(2) or "16.8", "--ipc", "0",
                "--cm-beta", f"{b:.5f}"], "pairs", None
    if m := re.fullmatch(r"(M1|E1)_En([\d.]+)", sample):
        w, b = w_beta(float(m.group(2)))
        return ["--energy", f"{w:.4f}", "--ipc", "1", "--ipc-multipole", m.group(1),
                "--ipc-min-open", "40", "--cm-beta", f"{b:.5f}"], "pairs", None
    if sample in GUN:
        return [], "neutrons", sample
    sys.exit(f"unknown sample '{sample}'")


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--run", required=True, help="B1 (pairs), B2 (beam n), B3 (beam γ), ...")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--target", required=True, choices=list(TARGETS))
    ap.add_argument("--det", required=True, choices=list(DETS))
    ap.add_argument("--njobs", type=int, default=10)
    ap.add_argument("--nevents", type=int, default=100_000)
    ap.add_argument("--flavour", default="workday")
    ap.add_argument("--memory", type=int, default=2048)
    ap.add_argument("--exe", required=True, help="frozen binary in bin/")
    ap.add_argument("--keep-root", type=int, default=2, help="the first N jobs also copy their ROOT file to EOS")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("extra", nargs=argparse.REMAINDER, help="-- extra simulation args")
    a = ap.parse_args()

    extra = a.extra[1:] if a.extra[:1] == ["--"] else a.extra
    name = f"{a.sample}_{a.target}_{a.det}"
    outdir = Path(EOS_BASE) / a.run / name
    jobdir = Path(JOB_BASE) / a.run / name
    exe = Path(a.exe).resolve()
    if not exe.is_file():
        sys.exit(f"ERROR: no binary {exe}")
    gen, prefix, beam = sample_args(a.sample)
    sim = target_args(a.target, beam) + gen + DETS[a.det] + extra

    seed0 = zlib.crc32(f"{a.run}/{name}".encode())
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
    if prefix == "neutrons":
        reduce = textwrap.dedent(f"""\
            mkdir -p "$D/parts_beam" "$D/parts_chan"
            ( set +u; source {LCG} > /dev/null
              python3 {ANALYSIS}/nfs_beam_reduce.py "${{W}}_t0.root" -o "$D/parts_beam/${{S}}_t0.npz"
              python3 {ANALYSIS}/bham_chan_reduce.py "${{W}}_t0.root" -o "$D/parts_chan/${{S}}_t0.npz" )
        """)
    else:
        reduce = textwrap.dedent(f"""\
            mkdir -p "$D/parts_nfs" "$D/parts_time"
            ( set +u; source {LCG} > /dev/null
              python3 {ANALYSIS}/nfs_reduce.py "${{W}}_t0.root" -o "$D/parts_nfs/${{S}}_t0.npz"
              python3 {ANALYSIS}/lnl_time.py "${{W}}_t0.root" -o "$D/parts_time/${{S}}_t0.npz" )
        """)
    run = textwrap.dedent(f"""\
        D=$(dirname "$OUT"); S=$(basename "$OUT")
        W="${{_CONDOR_SCRATCH_DIR:-/tmp/$USER}}/$S"
        "{exe}" -t 1 -n {a.nevents} -o "$W" -s "$SEED" {simargs} > "$OUT.log" 2>&1
        tail -4 "$OUT.log"
    """) + reduce + textwrap.dedent("""\
        if [ "$KEEP" = 1 ]; then cp "${W}_t0.root" "$D/"; fi
        rm -f "$W"*
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
    ] + [f"  {o}, {s}, {t}, {int(i < a.keep_root)}" for i, (o, s, t) in enumerate(jobs)] + [")"]
    sub.write_text("\n".join(lines) + "\n")
    info = (f"{a.run}/{name}\n{a.njobs} x {a.nevents} {prefix}\nmaster seed {seed0}\nexe {exe}\n"
            f"args: {' '.join(sim)}\nreduced in the job; ROOT kept for the first {a.keep_root} jobs\n")
    (jobdir / "README").write_text(info)
    (outdir / "RUN_INFO.txt").write_text(info)
    for i in range(4):     # condor_submit fails now and then from some nodes
        if os.system(f"condor_submit {sub}") == 0:
            return
    sys.exit("condor_submit failed")


if __name__ == "__main__":
    main()
