# MX17 simulation campaign — status & next steps

**Updated:** 2026-06-10 (late evening) · companion docs:
[PLAN_NEUTRON_CAMPAIGN.md](PLAN_NEUTRON_CAMPAIGN.md) ·
[docs/he3_self_shielding_note.md](docs/he3_self_shielding_note.md)

---

## Where things stand

### Jobs in flight (submitted 2026-06-10 ~23:40 CERN, 300 jobs, all queued healthy)

| Batch | Primaries | Output (EOS `mx17_geant_sim_results/`) | Purpose |
|---|---|---|---|
| Run B — sub-keV neutrons | 100 × 10M = 10⁹ | `neutrons_subkev/` | capture budget, wall backgrounds, singles, thermal-statistics check |
| Run B-full — full-range neutrons (1 meV–100 MeV) | 100 × 5M = 5×10⁸ | `neutrons_fullrange/` | flux-weighted picture across the board, MeV-region rates |
| Run A — X17+IPC pairs, new STEP geometry | 100 × 100k = 10M | `pairs_v2_step_target/` | new-geometry event pools, response JSON, MS budget |

Condor dirs (logs/submit files): `/afs/cern.ch/user/d/dneff/condor/mx17_{neutrons_subkev,neutrons_fullrange,pairs_v2}/`.
Repo on lxplus: `/afs/cern.ch/user/d/dneff/work/git/MX17_Full_Geant` (built clean, Geant4 11.2/LCG).

### New simulation capabilities (this session)

- **`--neutron <flux.root> <lambda2d.root>`** (event_type 2): EAR2 Ph3 evaluated
  flux + energy-dependent radial profile (`data/`); primary neutron's terminal
  interaction recorded in EventTree (`neutron_E_eV`, `capture_vol`,
  `capture_proc`, `cap_x/y/z`). Normalisation anchor: **7.31×10⁶ n/pulse < 1 keV**.
- **`--gamma-source <capture_lib.csv>`** (event_type 3): biased wall-background
  γ generator (IAEA-PGAA ²⁸Al/¹³C/²H cascade tables) fed by
  `scripts/make_capture_library.py` from a neutron run.
- `scripts/submit_neutrons.py`: Condor submission for both modes.

### Two Geant4 bugs found & fixed tonight (smoke tests on lxplus)

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
| ³He(n,p)t in gas | 73% | the dominant fate, as expected |
| escaped world | 24% | beam-profile tail outside 20 mm bore + window-top transmission |
| **Al vessel nCapture** | **7.6×10⁻³/n** | **~60× the analytic toy** — the ~5 mm on-axis dome + thick shoulder |
| LS/PVT H-capture | ~1.2% | 2.22 MeV γ born *inside* the calorimeter — background channel not in any toy |
| He-3 radiative (n,γ) | ~10⁻⁸ | ⇒ only ~10 events expected in all of run B: the He3-IPC channel **must** come from the explicit generator, never from neutron runs |

### The normalization issue (biggest open item — discuss with Alberto)

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

## Tomorrow

1. **Check jobs** (`condor_q`; expect overnight completion). Spot-check one
   file per batch: capture budget per volume/decade, no held jobs, file sizes.
2. **Capture budget vs Alberto's table & the analytic cap** — first physics
   deliverable of run B; directly tests the self-shielding correction.
3. **Run C**: `make_capture_library.py` on run-B output → 
   `submit_neutrons.py --gamma-source capture_lib.csv` (biased wall-background
   statistics).
4. **Run D**: single-particle KE scans (`--single`) for clean trigger-curve
   extraction P(trig | KE, θ) to refresh the response JSON.
5. **`analyze_neutrons.py`**: capture maps, singles rates per detector,
   conversion-pair spectra, trigger classification from runs B/B-full.
6. **Regenerate the response JSON + Highland budget from run A**
   (new geometry) and rerun the significance studies with it.
7. **Talk to Alberto**: thin-target vs self-shielding; what anchored the
   sensitivity; ROI choice (sub-keV vs MeV region); IPC coefficient
   provenance (2.1×10⁻³ — multipolarity? E0 contribution?).
8. If the MeV ROI wins: extend the pair generator with E_n-dependent
   kinematics and re-plan run A statistics accordingly.
9. Start `make_event_pools.py` (slimming) + the Python pile-up sampler
   (PLAN stages 3–4).
