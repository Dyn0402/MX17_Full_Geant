# Handoff

## B4C slab (July target) — updated 2026-10-20 (dylan-Yoga)

**Resume:** Done; used by nTof_x17/ntof_b4c_target. Nothing pending here.

**Goal:** simulate the 15 July B4C disc in the EAR2 beam to see if it could give pairs (it can't).

**Done:**
- `--slab-radius <mm>`, `--bias-slab <f>` (nCapture bias in the slab), and cascade entries in
  `--gamma-lines` (`E1+E2+...:w` emits every γ of the drawn entry from one vertex). Commit 1052a23.

**Gotchas / decisions:**
- G4 11.2 ParticleHP ¹⁰B(n,γ): γ lines sampled independently, 47 % of cascades sum above Q = 11.45 MeV.
  `/process/had/particle_hp/use_photo_evaporation true` (macro, before beamOn) conserves energy but its
  branching is statistical (~1/3 via a 2.12 MeV level not in EGAF). Use the EGAF cascade via --gamma-lines.
- ¹⁰B(n,α) is neutronInelastic, so --bias-slab only boosts (n,γ). The 478 keV ⁷Li line comes in ~56 % of (n,α), not 94 %.
- `-t >1` aborts ("gInterpreter not initialized") in the lxplus LCG build: use -t 1.
- lxplus build/worktree: `/afs/cern.ch/work/d/dneff/git/MX17_Full_Geant_b4c`; outputs `full_sim/b4c_july/`.

## Gas primary-ionization fix + capsule-Al capture bias — 2026-10-15 (dylan-MS-7C84)

**Resume:** ported and validated locally; build on lxplus in a new clone (the ear2 clone is in use by running jobs), then submit.

- Port of MX17_Geant e17c72c: `GasProperties` (ideal gas at 20 C, volume -> mass fractions, Garfield W/Fano), PAI in a
  `GasRegion` (DriftGas + AmpGas) via `G4EmConfigurator` with the region's own `G4ProductionCuts`, fluo/Auger/ignore-cut.
  Validated (Geant4 11.3.2, 20k cosmics, same seed vs 19a6a71): ArIso 1.698 mg/cm3 (ideal 1.699; was 1.826); gas Edep per
  (event, arm) with a gap-crossing muon, median per-event ratio new/old 0.906 (MX17_Geant: 0.908); muon step median
  47 -> 26 eV (PAI active, not a no-op); scintillators unchanged (0.98-1.02); CPU unchanged.
- `--bias-capsule <f>`: nCapture occurrence bias on `He3Cap_Al` only (weight 1/f in EventTree.weight). Validated at f = 20,
  200k EAR2 neutrons: Al captures 9861 raw / 575 weighted vs 579 analog; 3He weighted 89301 vs 89464; Al-capture events
  with an A/C gap track 444/384 vs 22/26 (x17); CPU +12 %; files x2. Capture positions shift slightly toward the nose
  (r < 8 mm: 0.55 vs 0.50), as expected from front-loading; weights restore rates.
- nTof_x17 side: `g4_digi` converts gas Edep with W = 26.0 eV and normalises ADC/e to 247 e- per track measured on the
  OLD sims (`run_digi.py:145`); both need updating for libraries from this build.


## Branch merges and build check — updated 2026-10-14 (dylan-MS-7C84)

**Resume:** Merged ear2-beam, ill, ill_ring, lnl and trigger_plastics into main (pushed, tip dbe53d0). The build is NOT verified: Geant4 is not installed on dylan-MS-7C84.

**Goal:** all of this repo's feature branches on main, and the code still building.

**Done:**
- Merges into main, one at a time, no conflicts. `ill`, `ill_ring` and `trigger_plastics` were already contained in `ear2-beam` and `lnl`.
- Local copies of those branches deleted.

**In progress / where it stopped:**
- `build/` is stale: its CMakeCache points at a different directory (`build2`), so `cmake --build build` fails with "Makefile: No such file".
- A fresh configure in a temp directory fails at `find_package(Geant4)`: Geant4 is not installed here.

**Next steps:**
1. On a machine with Geant4: `cmake -S . -B build-new && cmake --build build-new -j8`, then run the macros.
2. Once that builds, decide whether to delete the stale `build/` (not done here).

**Gotchas / decisions:**
- The merged branches all forked from 15bd0c2; `lnl` was the checked-out branch at the time.

**Key files & commands:**
- `CMakeLists.txt` — needs Geant4 on the CMake path

## Scattering-budget upgrade note — updated 2026-10-06 (dylan-MS-7C84)

**Resume:** Note is published (dylan-neff.web.cern.ch/notes/scattering-upgrade.html); next: replace the Highland-ratio scaling with a full Geant4 run of the 2 bar mylar cell.

**Goal:** Show colleagues how much e+e- scattering (and X17 opening-angle smearing) drops if the 500 bar capsule becomes the 2 bar mylar-wrap cell (x17_facility_search/vessel_design) and the MM drift window becomes two aluminised mylar foils. Appendix: ⁴He bag instead of the 24 cm air gap.

**Done:**
- Found the baseline study: `docs/angular_resolution/angular_resolution_note.{md,pdf}`.
- `scripts/scattering_geometry.py`: ray-traced, acceptance-averaged x/X0 per component (real capsule profile from DetectorConstruction.cc, 4 flat MM arms at 250 mm, hexagonal cell skin).
- `scripts/make_scattering_upgrade_figs.py`: matplotlib figs 1-3 + appendix (`docs/scattering_upgrade/*.png`, `results.json`).
- `scripts/make_scattering_upgrade_deck.py`: slide note (slidedoc), published with `add-note.py --slug scattering-upgrade --force --deploy`.
- Results: x/X0 1.84% -> 0.13% (14x); sigma68(opening angle) 15.6 -> 3.8 deg; peak height x0.22 -> x0.50; with 4He bag 2.5 deg.
- slidedoc (dylan-cern-site) gained left/right arrow slide navigation.

**In progress / where it stopped:** nothing half-done; the published note matches the scripts.

**Next steps:**
1. Full Geant4 simulation of the 2 bar cell geometry to replace the Highland-ratio scaling of the pairs_v2 response.
2. Model the upgrade acceptance (300 mm cell has a long vertex distribution; soft legs the old wall stopped now get through).
3. Decide the real window foil thicknesses (12 um each is an assumption) and whether the 4He bag is worth it.

**Gotchas / decisions:**
- Baseline psi tables are the pairs_v2 'first' estimator (direction at first MM hit, no vertex constraint) - used for both cases because the larger cell loses the target-centre vertex trick. Drift gas is included in every figure.
- Ratio scaling uses mean x/X0 per track; baseline reproduces the note's sigma68 = 15.5 deg (check passed).
- About 25% of accepted baseline tracks exit through dome/neck (2.36% X0 vs 1.36% barrel); 0.1% of upgrade tracks hit end caps (dropped).
- Cell rods/seam (~9% of azimuth) ignored; 7 um Al permeation layer would take sigma68 3.8 -> 3.9 deg.
- Plot label says "scattering", not "median scattering", because sigma68 is not a median.

**Key files & commands:**
- `python3 scripts/make_scattering_upgrade_figs.py` - budgets, figures, results.json
- `python3 scripts/make_scattering_upgrade_deck.py` - rebuild the note HTML
- `python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py docs/scattering_upgrade/scattering-upgrade.html --slug scattering-upgrade --force --deploy`
