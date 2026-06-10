# Handoff: Full MX17 Experiment Geant4 Simulation

## Context

This handoff is for a **new Claude session** that will build a new Geant4 repository
for the full MX17 detector simulation.  Three codebases are relevant — read all three
before writing any code.

| Codebase | Path | Purpose |
|----------|------|---------|
| Existing Geant4 (single-arm) | `/home/dylan/CLionProjects/MX17_Geant` | Mature single-detector-arm sim; materials, physics list, detector stack, analysis scripts all done |
| Python fast MC | `/home/dylan/PycharmProjects/nTof_x17/MX17_Simulation` | Acceptance / trigger / coincidence study; defines full 4-arm geometry and DAQ logic |
| X17 kinematics script | `/home/dylan/PycharmProjects/nTof_x17/e=e-_rel_angle_sim.py` | Generates X17 → e⁺e⁻ pairs with correct relativistic kinematics; opening angle, energy splitting |

---

## Physics Goal

The MX17 experiment searches for the X17 boson via the reaction **³He + n → ⁴He + X17**
followed by **X17 → e⁺e⁻**.  The detector surrounds a pressurised He-3 target on 4 sides
and reconstructs the e⁺e⁻ pair opening angle (~120°) and invariant mass.

The new simulation must reproduce the **full detector acceptance, angular resolution,
and calorimetric energy measurement** for signal pairs, so that the Python fast-MC
can sample from it rather than using a crude parameterisation.

---

## Simulation Strategy (agreed)

### Tier 1 — Full Geant4 (new repo): X17 signal only
- Fire **correlated e⁺e⁻ pairs** from the He-3 target centre with X17 kinematics
  (opening angle distribution from the kinematics script, proper energy splitting)
- Track both particles through the **4-arm detector stack** with full material budget
  and physics (EM + hadronic, same physics list as the single-arm sim)
- Output: **per-hit ROOT tree** (detector arm ID, sub-detector, local hit position u/v,
  time, edep) — one entry per hit, not per event
- Run large statistics (~10⁶–10⁷ pairs) once; store the hit tree permanently
- This gives correct correlated pair acceptance, angular resolution, and calorimetry
  including all material-budget effects

### Tier 2 — Python fast-MC (existing repo): backgrounds + analysis loop
- **IPC background pairs**: exponential opening-angle spectrum, generated cheaply in
  Python (straight-line tracks, geometric acceptance only)
- **Single-track random background**: isotropic, Poisson rate, generated in Python
- **X17 signal events**: **sampled from the Tier 1 ROOT hit tree** — pick a random event
  from the Geant4 output and inject its hits into the Python trigger window
- The coincidence trigger logic (plastic scint + LS1 coincidence on ≥2 arms),
  hit merging, angular reconstruction, and calorimetry reconstruction all live in Python
  and operate on the combined hit list (Geant4 X17 hits + Python background hits)
- This gives fast reruns (tweak rates, trigger thresholds, backgrounds) without
  rerunning the expensive Geant4 sim

**Why this split is correct:**
The pair acceptance cannot be factored into independent per-track efficiencies because
the azimuthal orientation of the pair matters (both tracks must hit *opposite* arms).
Only a 4-arm Geant4 run correctly captures this correlation.  Once captured in the hit
tree it can be sampled indefinitely.

---

## 4-Arm Geometry (from the Python simulation)

The He-3 target is at the origin.  Four identical detector arms are placed at
**+X, −X, +Y, −Y** (the beam is along Z; arms are in the transverse plane).

Default distances and sizes (from `run_example.py` — treat as starting point,
user may want to update):

| Component | Distance from target centre | Active area |
|-----------|----------------------------|-------------|
| MM drift gas front face | 22.0 cm | 38 × 34 cm |
| Plastic scintillator | ~28.0 cm | 48 × 48 cm |
| LS layer 1 front face | ~31.0 cm | 38 × 38 cm |

The arm ordering in the Python sim:
- Side 0: +X (detector normal points −X)
- Side 1: −X (detector normal points +X)
- Side 2: +Y (detector normal points −Y)
- Side 3: −Y (detector normal points +Y)

**Important:** the plastic scintillator is *larger* than the MM and LS to catch
tracks that scatter outward after the MM.  Model this correctly.

---

## Full Detector Stack per Arm (from single-arm sim)

Each arm is the **same material stack** as already built in the single-arm Geant4 sim.
Re-use these definitions exactly (copy the relevant DetectorConstruction code):

```
He-3 target (cylinder, axis ⊥ beam, already modelled)
  ↓  [air gap to arm front face]
GasWindow_Mylar   40 µm
GasWindow_Al       0.1 µm
DriftCathode_Kapton 50 µm
DriftCathode_Cu     9 µm
DriftGas           30 mm  ← primary MM tracking volume (100 µm step limit)
Micromesh          30 µm
AmpGas            150 µm
ResistivePaste    100 µm
PCB_Kapton         50 µm
PCB_Cu ×4          26 µm each
PCB_FR4 ×4        100 µm each
PCB_Rohacell        5 mm
PCB_AlFoil         50 µm
[air gap 20 mm]
ScintWall_BlackTape1  165 µm
PlasticScint           3 mm  ← trigger scintillator
ScintWall_BlackTape2  165 µm
ScintWall_AlFoil       50 µm
[air gap 20 mm]
LS_CFRP_1   1.5 mm  }
LiqScint_1  15 mm   }  ← calorimeter
LS_CFRP_2   1.5 mm  }
LiqScint_2  15 mm   }
LS_CFRP_3   1.5 mm  }
LiqScint_3  15 mm   }
LS_CFRP_4   1.5 mm  }
LiqScint_4  15 mm   }
LS_CFRP_5   1.5 mm  }
```

All material definitions (He-3 at 500 bar, CFRP, FR4, Rohacell 51, LAB, resistive
paste, PVC tape) are in `MX17_Geant/src/DetectorConstruction.cc` `DefineMaterials()`.
Copy that function verbatim.

The same physics list (`MX17_Geant/src/PhysicsList.cc`) should be used:
`G4EmStandardPhysics_option4` + `FTFP_BERT_HP` + step limiter.

---

## What to Score (output hit tree)

The output should be a **flat hit-level ROOT tree** (one entry per hit), NOT per-event
edep sums like the single-arm sim.  This lets the Python sim reconstruct tracks and
apply trigger logic.

Suggested branches per hit:

| Branch | Type | Description |
|--------|------|-------------|
| `eventID` | Int | G4 event number (links hits from same pair) |
| `trackID` | Int | G4 track ID (1=primary e⁻, 2=primary e⁺, >2=secondary) |
| `parentID` | Int | G4 parent track ID |
| `armID` | Int | 0=+X, 1=−X, 2=+Y, 3=−Y |
| `detType` | Char[32] | "DriftGas", "PlasticScint", "LiqScint_1", etc. |
| `particle` | Char[32] | "e-", "e+", "gamma", etc. |
| `u`, `v` | Double | Hit position in detector local frame [mm] |
| `w` | Double | Depth into sensitive volume [mm] |
| `edep` | Double | Energy deposited [eV] |
| `ke` | Double | Particle kinetic energy at step start [MeV] |
| `time` | Double | Global time [ns] |
| `gx`,`gy`,`gz` | Double | Global hit position [mm] |
| `px`,`py`,`pz` | Double | Particle momentum direction (unit vector) |

Also keep a per-event tree with the **generated pair kinematics** (both tracks' true
energy, direction, opening angle) for efficiency denominator calculations.

Score hits in: DriftGas, PlasticScint, LiqScint_1–4 (the four active volumes per arm).
Optionally score edep in AmpGas and the resistive paste (for MM gain studies).

---

## X17 Primary Generator

The kinematics script at `/home/dylan/PycharmProjects/nTof_x17/e=e-_rel_angle_sim.py`
generates e⁺e⁻ pairs.  Read it carefully before implementing the Geant4 generator.

Key physics to reproduce:
- X17 mass: ~17 MeV (check the kinematics script for exact value used)
- X17 is produced at rest in the He-3 target (or with a small boost from the
  ³He(n,p)T kinematics — check if the script accounts for this)
- The opening angle distribution between e⁺ and e⁻ is peaked near 120° but has
  a spread depending on the boost
- Energy is not equal-split — sample the correct energy asymmetry

Implementation: `G4VUserPrimaryGeneratorAction` with a `G4ParticleGun` fired twice
per event (once for e⁻, once for e⁺).  Store the generated kinematics in a
thread-local struct and write to the per-event tree.

Also implement a **single-particle mode** (fire one e⁻ or e⁺ at a chosen energy and
angle) for efficiency-vs-energy cross-checks against the single-arm sim.

---

## He-3 Target

The He-3 pressurised cylinder is already modelled in the single-arm sim:
- **Cylinder**: 3 cm diameter, 5 cm long, axis along Y (beam direction)
- 500 bar, 62.7 mg/cm³, pure ³He (G4Isotope)
- Outer walls: 0.5 mm Al + 1.2 mm CFRP
- The **gun fires from the centre** of the He-3 gas volume

Copy the He-3 cylinder geometry directly from
`MX17_Geant/src/DetectorConstruction.cc` (the `kFullExperiment` branch).

---

## Suggested Repo Structure

```
MX17_Geant_Full/
├── CMakeLists.txt
├── mx17_full_sim.cc          # main, CLI: -n events, -o outfile, -s seed, -t threads
├── include/
│   ├── SimConfig.hh          # geometry parameters (distances, sizes)
│   ├── DetectorConstruction.hh
│   ├── DetectorArm.hh        # single arm, placed 4 times with rotations
│   ├── PhysicsList.hh        # copy from single-arm sim
│   ├── ActionInitialization.hh
│   ├── X17PrimaryGenerator.hh
│   ├── HitData.hh            # per-hit data struct
│   ├── EventData.hh          # per-event (generated kinematics + accumulators)
│   ├── SteppingAction.hh
│   ├── EventAction.hh
│   └── RunAction.hh
├── src/
│   ├── DetectorConstruction.cc  # world + He-3 target + 4×DetectorArm
│   ├── DetectorArm.cc           # one arm geometry (reused 4 times)
│   ├── PhysicsList.cc           # copy from single-arm sim
│   ├── X17PrimaryGenerator.cc   # pair kinematics
│   ├── SteppingAction.cc        # fill HitData per step
│   ├── EventAction.cc
│   └── RunAction.cc             # write hit tree + event tree
└── scripts/
    ├── analyse_acceptance.py    # reads hit tree, computes pair acceptance vs angle
    └── sample_x17_hits.py       # helper for Python fast-MC to sample Geant4 hits
```

---

## Integration with Python Fast-MC

The Python sim (`MX17_Simulation/MX17_Simulator.py`) expects hits as objects with:
- detector_id, detector_type (string)
- position (u, v in local frame)
- time
- particle type
- pair_id (to link both particles of a pair)

Write `scripts/sample_x17_hits.py` as a bridge:
```python
def sample_x17_event(hit_tree, rng) -> list[Hit]:
    """Pick a random Geant4 event and return its hits in Python sim format."""
```

The Python sim then calls this instead of generating X17 pairs geometrically,
getting all the correct material-budget effects for free.

---

## Key Things to Get Right

1. **Local coordinate frame per arm**: each arm needs a local (u, v, w) frame
   where u and v are the two transverse directions and w is depth.  The Python
   sim uses this for spatial smearing and clustering.

2. **Arm rotation**: arms at ±X are rotated 90° around Z relative to arms at ±Y.
   Place them using G4RotationMatrix so the local u-axis is always horizontal
   (parallel to the beam axis Z) and v is vertical.  Confirm with the Python sim's
   coordinate convention.

3. **Transverse detector sizes**: plastic scint (48×48 cm) is larger than MM
   (38×34 cm) — this is intentional and must be reproduced.

4. **Step limit in gas**: 100 µm in DriftGas and AmpGas (same as single-arm sim).

5. **Thread safety**: use per-thread HitData accumulation (same pattern as single-arm
   sim's EventData).  Each thread writes its own output file; merge with `hadd`.

6. **Geometry check**: run with `/vis/viewer/flush` or write a sanity-check script
   that verifies all 4 arms are correctly placed and symmetric.

---

## Do NOT Reproduce from Scratch

These are already correct in the single-arm sim — copy verbatim:
- All material definitions in `DefineMaterials()`
- The physics list (`PhysicsList.cc`)
- The `SensitiveDetector` class (lightweight, scoring is in SteppingAction)
- The ROOT output pattern in `RunAction` (per-thread files, `hadd` to merge)
- The `ActionInitialization` pattern
- The `SimConfig` struct approach (just extend it with new geometry parameters)

---

## Questions to Clarify with Dylan Before Starting

1. **X17 mass and production kinematics**: does the generator account for the
   He-3 + n → He-4 + X17 boost, or is X17 produced at rest?  Read
   `e=e-_rel_angle_sim.py` and confirm.

2. **Detector distances**: the Python sim uses 22 cm MM face distance.  Is this
   fixed or should it be a config parameter?

3. **Number of LS layers per arm**: Python sim only mentions `liq_scint_1` as a
   trigger layer.  Should all 4 LAB layers be simulated per arm, or just the first?

4. **Beam axis**: is Z the beam axis (n-TOF beam direction)?  Confirm the arm
   placement is in the XY plane (transverse to beam).

5. **Magnetic field**: none mentioned in Python sim — confirm not needed.
