# MX17 simulation campaign — status & next steps

**Updated:** 2026-06-12 · companion docs:
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
