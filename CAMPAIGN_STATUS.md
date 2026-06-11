# MX17 simulation campaign — status & next steps

**Updated:** 2026-06-11 (morning) · companion docs:
[PLAN_NEUTRON_CAMPAIGN.md](PLAN_NEUTRON_CAMPAIGN.md) ·
[docs/he3_self_shielding_note.md](docs/he3_self_shielding_note.md)

---

## Current state (2026-06-11 ~11:15 CERN)

### All 300 original jobs completed overnight

| Batch | Files landed | Total size | Notes |
|---|---|---|---|
| Run B — sub-keV neutrons | **100 / 100** | 188 GB | Complete |
| Run B-full — full-range neutrons | **98 / 100** | 126 GB | jobs 076, 083 missing |
| Run A — X17+IPC pairs, new STEP geometry | **98 / 100** | 215 GB | jobs 030, 078 missing |

The 4 missing files all failed with **disk quota exceeded** on the personal EOS area — the
simulation ran to completion (5M / 100k events written) but the output ROOT file could
not be opened.  No partial files exist; re-running is clean.

---

## In progress right now

### 1. EOS migration (background copy running on lxplus)

Output has been relocated from the personal user area to the shared n_TOF experiment space:

```
OLD: /eos/user/d/dneff/mx17_geant_sim_results/{neutrons_subkev,neutrons_fullrange,pairs_v2_step_target}/
NEW: /eos/experiment/ntof/data/x17/full_sim/{neutrons_subkev,neutrons_fullrange,pairs_v2_step_target}/
```

`submit_neutrons.py` and `submit_pairs.py` defaults updated to the new path and pushed.
The copy is running as a `nohup` process on lxplus; log at
`/afs/cern.ch/user/d/dneff/eos_copy.log` (~530 GB total; may take an hour or more).

**After the copy completes:**
```bash
# verify counts
ls /eos/experiment/ntof/data/x17/full_sim/neutrons_subkev/ | wc -l   # expect 100
ls /eos/experiment/ntof/data/x17/full_sim/neutrons_fullrange/ | wc -l # expect 98
ls /eos/experiment/ntof/data/x17/full_sim/pairs_v2_step_target/ | wc -l # expect 98
# then remove originals to free ~530 GB of personal quota
rm -r /eos/user/d/dneff/mx17_geant_sim_results/{neutrons_subkev,neutrons_fullrange,pairs_v2_step_target}
```

### 2. Four retry jobs in Condor queue (writing directly to ntof EOS)

| Condor cluster | Jobs | Batch | Seeds |
|---|---|---|---|
| 16248995 | 076, 083 | neutrons_fullrange | 1850501473, 1420052173 |
| 16248996 | 030, 078 | pairs_v2_step_target | 1170252924, 1239854304 |

Submit files: `/afs/cern.ch/user/d/dneff/condor/mx17_{neutrons_fullrange,pairs_v2}_retry/`
These write directly to the new ntof path, so nothing more to do when they finish.
Monitor with `condor_q dneff`.

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

## Next steps (once copy + retries finish)

1. **Verify & clean up**: check copy log, confirm file counts on ntof EOS, then
   `rm -r` the originals from `/eos/user/d/dneff/mx17_geant_sim_results/`.

2. **Capture budget** (`analyze_neutrons.py`) — first physics deliverable of run B:
   - Extract `capture_vol` + `capture_proc` + `neutron_E_eV` from run B sub-keV output.
   - Bin by volume (He3 gas, Al vessel, LS/PVT, CFRP, world escape) and energy decade.
   - Compare Al capture rate to analytic toy (~1.2×10⁻⁴/n expected) and smoke-test value (7.6×10⁻³/n).
   - Cross-check He3 (n,p)t rate vs Alberto's self-shielding-corrected table.
   - Flag the LS/PVT H-capture channel (2.22 MeV γ inside calorimeter) — quantify per detector.

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
