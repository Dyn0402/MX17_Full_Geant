# MX17 simulation campaign — status & next steps

**Updated:** 2026-07-23 · companion docs:

---

## ⚠ 2026-07-23: capsule flipped to NOSE-FIRST — all results below are stale

The mounting audit confirmed the He-3 capsule sits **tip into the beam**; the
sim had it valve-first from 2026-06-10 onward (inherited `rotateX(-90°)`,
fixed in `3d97437`). Al-vessel captures drop **7.81e-3 → 4.76e-3 per neutron**
(×0.61) and move from the valve stem to the nose, and the gas capture-depth
profile flips with it — so the thermal-gate background, the Al-pair/leg rates,
and the pair-vertex library all change. Every number on this page and in
`docs/report/*`, `docs/al_gamma_yield_check/`, and the thermal/MeV notes was
produced valve-first.

Re-running (nose-first, new `*_2cm_nose` EOS dirs; the valve-first data is
kept for comparison):

| dataset | jobs × events | status |
|---|---|---|
| `neutrons_thermal_trig_2cm_nose` | 100 × 10M, 1 meV–2 eV | submitted 07-23 |
| `neutrons_epi_trig_2cm_nose` | 50 × 10M, 2 eV–100 keV | submitted 07-23 |
| `neutrons_thermal_bias1e5_2cm_nose` | 25 × 4M, ×1e5 ³He(n,γ) bias | submitted 07-23 |
| `pairs_thermal_trig_2cm_nose` | 100 × 100k | waits on the new gas-vertex library |
| `gsrc_mechanism_nose` | gamma-source, leg mechanism | waits on the new capture library |

Not re-run (closed cross-checks): `--no-al` and the γ-cut scans — 100 µm is
settled. The conv-pair study reuses `neutrons_thermal_trig_2cm_nose`, since
the current build writes `ConvPairTree` in every run.

---

## Headline result (2026-07-19): thermal-gate (>1 ms) trigger optimized — final geometry

Three campaigns with the FINAL surveyed geometry (flipped stack, STEP LS
vessels, pinwheel MMs; commits through `808ffe0`), all validated 250/250 on
`/eos/experiment/ntof/data/x17/full_sim/`:
`neutrons_thermal_trig` (10⁹, 1 meV–2 eV), `neutrons_epi_trig` (5×10⁸,
2 eV–100 keV), `pairs_thermal_trig` (10⁷ X17+IPC 50/50, vertices from the
measured thermal self-shielding profile via the new `--pair-vertex-lib`;
median capture depth 14 mm into the gas, 95 % in the first 25 mm — the pair
source sits ~3.7 cm upstream of the target centre).
Gate: >1 ms at 19.5 m ⇔ E_n < 1.99 eV, **4.284×10⁶ n/pulse in-gate**.
MIP calibration (20k µ⁻): **1 MIP = 458 keV (SiPM bar) / 4.33 MeV (plastic)**.
Analysis: `scripts/analyze_trigger_thermal.py` →
`analysis/trigger_thermal/` (lxplus): `trigger_scan.json`, spectra + ROC.

Findings (a = SiPM-bar, b = plastic-bar thresholds in MIP, per-channel):
- **Epithermal spill-in is negligible** (delayed H-capture leakage into the
  gate: 0 pair-tags in 5×10⁸ n; singles ~8/pulse vs 170 thermal). The gate
  background is genuinely thermal-arrival captures (Al dome 6.5×10⁻³/n,
  H-captures in LS 2.1×10⁻³/n and plastics ~1.4×10⁻³/n at thermal).
- Requiring the plastic in **both** legs costs ~5×: plastics subtend far
  less solid angle than the wall (30 cm vs 50 cm in v, 10 cm further back).
  ε(X17) ceiling: 3.1 % (2 full legs) vs 10.8 % (2 SiPM + ≥1 confirm) vs
  15.5 % (SiPM-only, but 58 bg pair-tags/pulse).
- Plastic spectrum physics: signal legs punch through 2.5 cm (crossing e±
  deposit ~1 MIP; stopped legs up to ~3 MIP); thermal-bg plastic hits die
  above ~1.7 MIP (7.72 MeV Al Compton edge), H-capture γ stays < 0.5 MIP.
  SiPM bg shows no MIP peak (γ Comptons) — signal peaks cleanly at 1 MIP.
- **Recommended menu: 2 SiPM legs ≥ 0.5 MIP + ≥1 plastic confirm at
  1.0–1.2 MIP (4.3–5.2 MeV, ≈ the toy study's 5–6 MeV per-leg cut):**
  ε(X17) = 7.0→4.0 %, ε(IPC, both legs >5 MeV) = 5.3→3.2 %, background
  1.8→0.8 pair-tags/pulse (0.4→0.2 Hz at 1.93×10⁴ pulses/day) — DAQ-trivial.
  b ≥ 1.5 MIP kills the crossing-leg signal (ε → 1.4 %); avoid.
- Full 2-leg coincidence as the high-purity alternative: (0.5, 0.5–0.8) MIP
  → ε 2.1–1.3 %, bg 0.17–0.03/pulse.
- Context: in-gate ³He(n,γ) sits on the self-shielding ceiling
  (~1×10⁻⁴ IPC/pulse), so the >1 ms gate is background-characterization
  territory; the same menu should be re-scanned for the µs-scale MeV window
  (rates ×10⁴ higher there — separate study needed).
[PLAN_NEUTRON_CAMPAIGN.md](PLAN_NEUTRON_CAMPAIGN.md) ·
[docs/he3_self_shielding_note.md](docs/he3_self_shielding_note.md) ·
[docs/report/thermal_note.pdf](docs/report/thermal_note.pdf) ·
[docs/report/mev_note.pdf](docs/report/mev_note.pdf)

---

## Headline result (2026-06-12): MeV region confirmed feasible — 22.5 X17/day in 0.2–2 MeV

Run B-full analysed (100 validated jobs, 5×10⁸ neutrons, 1 meV–100 MeV,
`scripts/analyze_mev_captures.py` + `make_mev_report_figures.py`): **714
direct ³He(n,γ) events** (95 % above 100 keV, median 1.3 MeV) give
**32.3 capt/pulse → 32.7 X17 produced/day** over the full range, of which
**22.5 ± 1.0/day in the 0.2–2 MeV window** — ×300 the sub-keV yield and
within 10 % of the commented `results_3He` row. Where thin-target is valid
the table is *confirmed* (per-neutron G4/table = 1.2–1.7); below 1 keV the
data sit on the opacity ceiling, reproducing the thermal result. Effective
σ_nγ/σ_np per decade tracks ENDF over six orders of magnitude. Window sits
1.0–3.2 µs after the γ flash at 19.5 m, sharing the gate with 3×10⁵ (n,p)t
and 9×10⁴ scintillator H-captures per pulse. Production is settled —
remaining work is acceptance (pairs run) and trigger/pile-up design.
Write-up: [docs/report/mev_note.pdf](docs/report/mev_note.pdf) (8 pp);
rates: `analysis/mev/mev_rates.json`.
**Normalisation note:** fullrange n/pulse = 2.263×10⁷ (flux-file integral;
do not reuse the 7.31×10⁶ sub-keV anchor, and Alberto's 3.29×10⁷ is 31 %
hotter — compare per neutron). **Physics note:** never apply the thermal
σ_nγ/σ_np = 10⁻⁸ above ~1 keV; use direct counts.

## Headline result (2026-06-11, FINAL — full 10⁹ statistics): thermal statistics confirmed dead

Run B complete (100 validated jobs, 10⁹ neutrons,
`scripts/analyze_thermal_captures.py`): the sub-keV ³He(n,γ) rate sits at
the self-shielding ceiling — **15 direct (n,γ) events** observed, all
eleven below 100 eV absorbed within mm of the gas entrance face (textbook
self-shielding). The direct counts also resolve the ENDF energy dependence:
ratio (2.0±0.5)×10⁻⁸ overall, thermal events at 1.0×10⁻⁸, eV–keV events at
the elevated ENDF σ_nγ/σ_np. Integrated: **(1.1–2.3)×10⁻⁴ IPC/pulse vs the
table's 1.21×10⁻² — factor ×50–100** (direct ×53, thermal-ratio floor ×106).
Sub-keV-anchored sensitivity scales 3.0σ → ~0.3–0.4σ. The MeV region
(thin-target valid, ~98 % of the rate) is the likely ROI.
Full write-up: [docs/report/thermal_note.pdf](docs/report/thermal_note.pdf)
(12-page internal note, final statistics, incl. ENDF cross-section figure
and measured absorption-position maps).

## Quota corruption post-mortem (worse than first thought)

The overnight quota exhaustion didn't just kill 4 jobs — every job that hit
the full personal-EOS quota **mid-write** produced a GB-sized ROOT file with
no key directory (`TBranch::WriteBasketImpl: WriteBuffer failed` in .err,
unreadable by uproot/ROOT). Validation (`EventTree.num_entries` check on all
files): **60/100 subkev, 14/100 fullrange (+2 pending first retries),
4/100 pairs corrupt**. All 78 resubmitted ~12:30 (clusters
16249004/5/6, original seeds, writing directly to ntof EOS).
**Lesson: validate every output file's tree readability after any campaign;
`ls`-level checks and byte counts are not enough.**

## Current state — campaign COMPLETE

All 300 files validated (tree readable, exact entry counts):
**100/100 each** of `neutrons_subkev/` (10⁹), `neutrons_fullrange/` (5×10⁸),
`pairs_v2_step_target/` (10⁷) at `/eos/experiment/ntof/data/x17/full_sim/`.
Originals deleted from personal EOS (~530 GB freed; usage 1.1 TB → 585 GB).
Scan results: `thermal_captures_subkev_full.npz/json` (repo on lxplus).

---

## Infrastructure summary

- **EOS output**: `/eos/experiment/ntof/data/x17/full_sim/`
- **Condor dirs** (logs/submit files): `/afs/cern.ch/user/d/dneff/condor/mx17_{neutrons_subkev,neutrons_fullrange,pairs_v2}/`
- **Repo on lxplus**: `/afs/cern.ch/work/d/dneff/git/MX17_Full_Geant` (built clean, Geant4 11.2/LCG)

---

## Physics context

### New simulation capabilities

- **`--neutron <flux.root> <lambda2d.root>`** (event_type 2): EAR2 Ph3 evaluated
  flux + energy-dependent radial profile (`data/`); primary neutron's terminal
  interaction recorded in EventTree (`neutron_E_eV`, `capture_vol`,
  `capture_proc`, `cap_x/y/z`). Normalisation anchor: **7.31×10⁶ n/pulse < 1 keV**.
- **`--gamma-source <capture_lib.csv>`** (event_type 3): biased wall-background
  γ generator (IAEA-PGAA ²⁸Al/¹³C/²H cascade tables) fed by
  `scripts/make_capture_library.py` from a neutron run.

### Two Geant4 bugs found & fixed (smoke tests on lxplus, 2026-06-10)

1. **³He(n,p)t is `neutronInelastic`, not `nCapture`** in the HP package —
   the capture hook initially recorded 0% of gas absorptions. Now records any
   terminal hadronic process, with the process name in `capture_proc`.
2. **`G4NeutronTrackingCut` removed from the PhysicsList** — its default 10 µs
   tracking limit killed slow neutrons mid-flight (0.4 eV neutron needs ~20 µs
   for 20 cm). Half of all sub-keV beam neutrons died as `nKiller` before
   reaching anything. *Any earlier results involving neutron transport carried
   this distortion.*

### Smoke-test physics (20k sub-keV neutrons, post-fix)

| Channel | Fraction | Comment |
|---|---|---|
| ³He(n,p)t in gas | 73% | dominant fate, as expected |
| escaped world | 24% | beam-profile tail outside 20 mm bore + window-top transmission |
| **Al vessel nCapture** | **7.6×10⁻³/n** | **~60× the analytic toy** — 5 mm on-axis dome + thick shoulder |
| LS/PVT H-capture | ~1.2% | 2.22 MeV γ born *inside* the calorimeter — not in any toy |
| He-3 radiative (n,γ) | ~10⁻⁸ | ~10 events expected in all of run B; He3-IPC must come from explicit generator |

### The normalisation issue (biggest open item — discuss with Alberto)

`calculation_tables/results_3He` computes radiative captures **thin-target**
(areal × σ_nγ), but below ~1 keV the gas is opaque to (n,p) (optical depth ~150
at 25 meV), capping radiative captures at σ_nγ/σ_np = 1.01×10⁻⁸ per incident
neutron. Flux-weighted sub-keV: table 1.21×10⁻² IPC/pulse → corrected
**1.9×10⁻⁴** (×62). If the ROI stays sub-keV, significance projections scale by
~1/8. The >100 keV rows are thin-target-valid and carry ~98% of the physical
rate — the likely resolution is a **MeV-region ROI** (note the commented
0.2–2 MeV row in the table: 5.2×10⁻² IPC/pulse), which would require
E_n-dependent generator kinematics (E* ≈ 20.58 MeV + ¾E_n + CM boost).
Full derivation: [docs/he3_self_shielding_note.md](docs/he3_self_shielding_note.md).

### Python fast-MC side (already in place, `nTof_x17/MX17_Simulation`)

Geant4-response smearing (`geant4_response.py` + old-geometry JSON), Asimov +
profiled significance machinery (`run_significance_study.py`), background toy
(`run_ipc_background_study.py`), Config A vs B report
(`results/significance/Config_A_vs_B_report.pdf`: A preferred, 3.0σ → 2.3σ
under data-driven background fit — all pre-self-shielding-correction).
Serial/parallel detector-plane inconsistency in `MX17_Simulator.py` fixed
(was ~40% optimistic MM acceptance in parallel mode).

---

## Next steps (once retries finish)

1. **Validate all 300 files** (tree readability + entry counts), then rerun
   `analyze_thermal_captures.py` on the full 10⁹ and refresh the numbers +
   figures in `docs/report/thermal_note.tex` (drop the PRELIMINARY banner).

2. **Capture budget done in the thermal scan** — Al 8.0×10⁻³/n (~60× toy),
   LS H-capture 1.2×10⁻² /n (largest non-gas channel). Remaining: per-detector
   breakdown and the run B-full MeV-region budget for the ROI decision.

3. **Run C** (biased wall-background statistics):
   ```bash
   python3 scripts/make_capture_library.py \
       /eos/experiment/ntof/data/x17/full_sim/neutrons_subkev/*.root -o capture_lib.csv
   python3 scripts/submit_neutrons.py --gamma-source capture_lib.csv \
       --outdir /eos/experiment/ntof/data/x17/full_sim/gammas_wall_bg \
       --nevents 1000000 --njobs 100
   ```

4. **Run D** — single-particle KE scans (`--single e-`) for trigger-curve extraction
   P(trig | KE, θ) to refresh the response JSON with the new STEP geometry.

5. **Regenerate response JSON + Highland MS budget from run A** (new geometry) and
   rerun significance studies.

6. **Talk to Alberto**: thin-target vs self-shielding; ROI choice (sub-keV vs MeV);
   IPC coefficient provenance (2.1×10⁻³ — multipolarity? E0 contribution?).

7. **If MeV ROI wins**: extend pair generator with E_n-dependent kinematics
   (E* ≈ 20.58 MeV + ¾E_n + CM boost) and re-plan run A statistics accordingly.

8. **`make_event_pools.py`** (Stage 3) + Python pile-up sampler (Stage 4).
