# MX17 Full Geant4 Simulation
## n_TOF X17 experiment — Dylan Neff

Full 4-arm Geant4 simulation of the MX17 detector for the reaction
**³He + n → ⁴He + X17**, **X17 → e⁺e⁻**.

The simulation fires correlated e⁺/e⁻ pairs from the He-3 target centre with
correct X17 kinematics (Lorentz-boosted decay) and tracks both particles
through the complete 4-arm detector stack. Output is a **flat hit-level ROOT
tree** (one entry per step in a sensitive volume) intended to be sampled by the
Python fast-MC (`nTof_x17/MX17_Simulation`) in place of geometric acceptance
approximations.

---

## Geometry overview

Beam axis = **+Y** (floor → ceiling). Four identical detector arms are placed
at **±X** (arms 0, 1) and **±Z** (arms 2, 3) in the transverse plane.
The He-3 target is at the **origin**, surrounded symmetrically by all four arms.

```
                  Arm 2 (+Z)
                     |
        Arm 1 (−X) ──●── Arm 0 (+X)       ● = He-3 target at origin
                     |
                  Arm 3 (−Z)

  Beam along +Y (into page)
```

Geometry plots (top-down and 3D) can be regenerated at any time:
```bash
python scripts/plot_geometry.py          # saves PNG files to scripts/
python scripts/plot_geometry.py --interactive   # opens live pyvista window
```

---

## He-3 target

Pressurised cylinder with axis along **Y** (beam direction).

| Component         | Material                   | Radius / thickness  |
|-------------------|----------------------------|---------------------|
| He-3 gas          | ³He at 300 bar, 37.6 mg/cm³ | r = 1.5 cm         |
| Al wall           | Aluminium                  | 0.5 mm              |
| CFRP outer shell  | C-fibre/epoxy, 1.55 g/cm³  | 0.9 mm              |

Cylinder length: **8 cm** along Y. Gun fires from the origin (He-3 gas centre).

---

## Detector stack per arm

The front face of the MM window is at `mm_distance` (default **22 cm**) from
the origin. All depths below are measured from that front face.

### Micromegas

| Layer                  | Material         | Thickness | G4 volume name        |
|------------------------|------------------|-----------|-----------------------|
| Gas window             | Mylar (PET)      | 40 µm     | `GasWindow_Mylar`     |
| Gas window coat        | Aluminium        | 0.1 µm    | `GasWindow_Al`        |
| Drift cathode          | Kapton           | 50 µm     | `DriftCathode_Kapton` |
| Drift cathode Cu       | Copper           | 9 µm      | `DriftCathode_Cu`     |
| **Drift gas** ★        | **active gas**   | **30 mm** | **`DriftGas`**        |
| Micromesh              | Stainless steel  | 30 µm     | `Micromesh`           |
| **Amp gas** ★          | **active gas**   | **150 µm**| **`AmpGas`**          |
| Resistive paste        | C/acrylic, 1.4 g/cm³ | 100 µm | `ResistivePaste`    |

Active area: **38 × 34 cm** (u × v), step limit 100 µm in gas volumes.

### PCB stack

| Layer        | Material              | Thickness | G4 volume name  |
|--------------|-----------------------|-----------|-----------------|
| Kapton       | Kapton                | 50 µm     | `PCB_Kapton`    |
| Cu × 4       | Copper                | 26 µm ea. | `PCB_Cu_1–4`    |
| FR4 × 4      | Epoxy-glass, 1.85 g/cm³| 100 µm ea.| `PCB_FR4_1–4` |
| Rohacell 51  | PMI foam, 0.052 g/cm³ | 5 mm      | `PCB_Rohacell`  |
| Al foil      | Aluminium             | 50 µm     | `PCB_AlFoil`    |

### Air gap 1 — 20 mm

### Trigger scintillator wall

| Layer           | Material              | Thickness | G4 volume name          |
|-----------------|-----------------------|-----------|-------------------------|
| Black mylar tape| Mylar (PET), light-tight | 200 µm | `ScintWall_BlackTape1`  |
| **Trigger scint ★** | **PVT plastic**   | **3 mm**  | **`PlasticScint`**      |
| Black mylar tape| Mylar (PET)           | 200 µm    | `ScintWall_BlackTape2`  |
| Al foil         | Aluminium             | 50 µm     | `ScintWall_AlFoil`      |

Active area: **48 × 48 cm** (u × v).

### Air gap 2 — 20 mm

### Liquid scintillator stack

Each LAB layer is preceded by an inner CFRP liner and Al reflector (inside the
structural CFRP walls). The full cell sequence is:

| Layer             | Material            | Thickness | G4 volume name      |
|-------------------|---------------------|-----------|---------------------|
| CFRP wall         | C-fibre/epoxy       | 2 mm      | `LS_CFRP_1`         |
| Inner CFRP liner  | C-fibre/epoxy       | 600 µm    | `LS_InnerCFRP_1`    |
| Al liner          | Aluminium           | 40 µm     | `LS_Al_1`           |
| **LS 1** ★        | **LAB, 0.86 g/cm³** | **2 cm**  | **`LiqScint_1`**    |
| CFRP wall         | C-fibre/epoxy       | 2 mm      | `LS_CFRP_2`         |
| Inner CFRP liner  | C-fibre/epoxy       | 600 µm    | `LS_InnerCFRP_2`    |
| Al liner          | Aluminium           | 40 µm     | `LS_Al_2`           |
| **LS 2** ★        | **LAB**             | **2 cm**  | **`LiqScint_2`**    |
| CFRP wall         | C-fibre/epoxy       | 2 mm      | `LS_CFRP_3`         |

Active area: **45 × 45 cm** (u × v).  Total stack depth: **~47.3 mm**.
LAB definition: C₁₈H₃₀ (JUNO recipe), density 0.86 g/cm³.

### Air gap 3 — 10 mm

### Back plastic scintillators

Two **20 × 30 cm × 2 cm** PVT bars per arm, placed side-by-side in the
horizontal transverse direction (u), giving **~40 cm (u) × 30 cm (Y/beam)**
combined coverage. Each bar is individually wrapped: Al foil directly on the
scintillator surface, then black mylar tape on the outside. ~3 mm gap between
wrapped bars.

| Component           | Material                 | Thickness   | G4 volume name     |
|---------------------|--------------------------|-------------|--------------------|
| Black mylar tape    | Mylar (PET), light-tight | 200 µm      | `BackScintTapeL/R` |
| Al foil             | Aluminium (reflector)    | 20 µm       | `BackScintAlL/R`   |
| **Left bar** ★      | **PVT plastic**          | **20×30×2 cm** | **`BackScintL`** |
| **Right bar** ★     | **PVT plastic**          | **20×30×2 cm** | **`BackScintR`** |

Geant4 nesting: tape envelope → Al envelope → PVT scint (all concentric).

★ = scored sensitive volume (hits recorded in output tree)

---

## Local coordinate frame per arm

Each arm has a local (u, v, w) frame used for hit positions in the output:

| Axis | Direction in world | Description |
|------|--------------------|-------------|
| u    | horizontal transverse (±Z for ±X arms, ±X for ±Z arms) | across the detector face |
| v    | +Y (beam axis)     | along the beam direction |
| w    | outward radial     | depth into the detector (0 = MM front face) |

Arm IDs: **0 = +X**, **1 = −X**, **2 = +Z**, **3 = −Z**.

---

## Physics list

| Component | Detail |
|-----------|--------|
| EM | `G4EmStandardPhysics_option4` — Livermore models below 100 keV, accurate Bremsstrahlung, Auger, delta rays |
| Hadronic | `FTFP_BERT_HP` — high-precision neutron data below 20 MeV |
| Decay | `G4DecayPhysics` + `G4RadioactiveDecayPhysics` |
| Step limits | 100 µm in Micromegas gas; 1 mm in He-3 gas |
| Production cuts | 10 µm for e±, 100 µm for γ/p |

---

## X17 primary generator

**Pair mode** (default): fires one e⁻ and one e⁺ per event from the origin.

- X17 mass: **16.8 MeV** (configurable via `--mass`)
- X17 lab energy: **20.58 MeV** (configurable via `--energy`), boosted along +Y
- Isotropic decay in X17 rest frame → Lorentz boost → correlated lab momenta
- Opening angle peaked near 120°; asymmetric energy split

**Single-particle mode** (`--single`): fires one particle at a fixed energy and
angle, useful for efficiency cross-checks against the single-arm simulation.

---

## Building on lxplus

```bash
ssh lxplus.cern.ch
git clone <repo> MX17_Full_Geant
cd MX17_Full_Geant

# Load Geant4 + ROOT from CVMFS (run once per login shell)
source scripts/setup_lxplus.sh

# Build (Release mode, up to 8 cores)
bash scripts/build.sh

# Optional clean rebuild
bash scripts/build.sh clean
```

The setup script automatically selects Geant4 11.2 (falling back to 11.1) and
ROOT via LCG on AlmaLinux 9. If a CVMFS node is not available the script will
print a clear error.

---

## Running

```bash
# X17 pairs, 10000 events, ArCF4 gas, 4 threads
build/mx17_full_sim -n 10000 -g ArCF4 -t 4 -o /eos/.../x17_run1

# Single-particle mode: 8 MeV e-, 90° from beam axis, toward arm 0 (+X)
build/mx17_full_sim -n 10000 --single e- 8 90 0 -o /tmp/singletrack_test

# Verbose (prints hit/event summary every 1000 events)
build/mx17_full_sim -n 1000 -v -o /tmp/test
```

### Command-line options

| Flag | Default | Description |
|------|---------|-------------|
| `-n <N>` | 10000 | Number of events |
| `-o <path>` | `x17_output` | Output file base (no extension) |
| `-s <seed>` | time-based | Random seed |
| `-t <N>` | 1 | MT threads (each writes its own file) |
| `-g <gas>` | `ArCF4` | MM drift gas: `ArCF4`, `ArIso`, `HeEth`, `ArCO2`, etc. |
| `-v` | off | Verbose event printout |
| `--mass <MeV>` | 16.8 | X17 mass |
| `--energy <MeV>` | 20.58 | X17 lab energy |
| `--dist <cm>` | 22.0 | Arm distance (MM front face to origin) |
| `--single <name> <E> <θ> <φ>` | off | Single-particle mode: `e-`, `e+`, etc. |

---

## Output

Two ROOT TTrees per file (`USE_ROOT` CMake flag, auto-detected).
In MT mode each thread writes `<outfile>_t<N>.root`; merge with `hadd`.

### HitTree — one entry per step in a sensitive volume

| Branch | Type | Units | Description |
|--------|------|-------|-------------|
| `eventID` | Int | — | Links hits from the same pair |
| `trackID` | Int | — | 1 = e⁻, 2 = e⁺, >2 = secondary |
| `parentID` | Int | — | Parent track ID |
| `armID` | Int | — | 0=+X, 1=−X, 2=+Z, 3=−Z |
| `detType` | Char[32] | — | `DriftGas`, `AmpGas`, `PlasticScint`, `LiqScint_1`, `LiqScint_2`, `BackScintL`, `BackScintR` |
| `particle` | Char[32] | — | `e-`, `e+`, `gamma`, … |
| `u`, `v`, `w` | Double | mm | Local hit position (u=transverse, v=beam, w=depth) |
| `edep` | Double | eV | Energy deposited in this step |
| `ke` | Double | MeV | Particle kinetic energy at step start |
| `time` | Double | ns | Global time |
| `gx`, `gy`, `gz` | Double | mm | Global hit position |
| `px`, `py`, `pz` | Double | — | Momentum unit vector |

### EventTree — one entry per event (generated pair truth)

| Branch | Type | Units | Description |
|--------|------|-------|-------------|
| `eventID` | Int | — | Event number |
| `vtx_x`, `vtx_y`, `vtx_z` | Double | mm | Production vertex in He-3 gas |
| `em_ke`, `ep_ke` | Double | MeV | e⁻ / e⁺ kinetic energy |
| `em_px/py/pz` | Double | — | e⁻ momentum unit vector |
| `ep_px/py/pz` | Double | — | e⁺ momentum unit vector |
| `openingAngle` | Double | deg | True opening angle |

The vertex is sampled uniformly within the He-3 cylinder (radius 1.5 cm, ±4 cm along Y).
Momentum unit vector + kinetic energy fully specify the 4-momentum of each lepton.
In single-particle mode the vertex is always (0, 0, 0).

### Merging MT output

```bash
hadd x17_merged.root x17_output_t*.root
```

### Reading in Python

```python
import uproot, numpy as np

with uproot.open("x17_merged.root") as f:
    hits = f["HitTree"].arrays(library="pd")
    evts = f["EventTree"].arrays(library="pd")

# Per-arm hit counts
print(hits.groupby("armID")["edep"].count())

# Events with hits in both LS layers on at least one arm
ls_arms = hits[hits.detType.isin(["LiqScint_1","LiqScint_2"])].groupby("eventID")["armID"].nunique()
print("Events with ≥2 LS-hit arms:", (ls_arms >= 2).sum())
```

---

## HTCondor submission (lxplus)

```bash
# Large X17 run — 1M events, 4 threads per job
python3 scripts/submit_condor_full.py \
    --outdir /eos/experiment/ntof/data/x17/geant_x17_hits \
    --jobdir /afs/cern.ch/user/d/dneff/condor/mx17_full \
    --nevents 100000 --threads 4

condor_q
```

After all jobs finish, merge per-thread files:

```bash
hadd /eos/.../geant_x17_hits/x17_merged.root \
     /eos/.../geant_x17_hits/x17_output_*.root
```

---

## Integration with Python fast-MC

The Python fast-MC (`nTof_x17/MX17_Simulation`) can sample X17 signal events
directly from the Geant4 hit tree instead of using geometric acceptance
approximations. See `scripts/sample_x17_hits.py` for the bridge helper:

```python
from scripts.sample_x17_hits import sample_x17_event

hits = sample_x17_event(hit_tree, rng)   # returns list of Hit objects
```

---

## File structure

```
MX17_Full_Geant/
├── CMakeLists.txt
├── mx17_full_sim.cc           Main: CLI parsing, RunManager setup
├── include/
│   ├── SimConfig.hh           All geometry + run parameters
│   ├── DetectorConstruction.hh
│   ├── HitData.hh             Per-hit data struct
│   ├── EventData.hh           Per-event struct (PairKinematics + hit vector)
│   ├── X17PrimaryGenerator.hh Pair kinematics / single-particle mode
│   ├── SteppingAction.hh      Fills HitData per step in scored volumes
│   ├── EventAction.hh
│   ├── RunAction.hh           Writes HitTree + EventTree
│   ├── ActionInitialization.hh
│   ├── PhysicsList.hh
│   └── SensitiveDetector.hh
├── src/
│   ├── DetectorConstruction.cc  He-3 target + 4 arms; all materials
│   ├── X17PrimaryGenerator.cc
│   ├── SteppingAction.cc
│   ├── EventAction.cc
│   ├── RunAction.cc
│   ├── ActionInitialization.cc
│   ├── PhysicsList.cc
│   └── SensitiveDetector.cc
├── macros/
│   └── run_default.mac
└── scripts/
    ├── setup_lxplus.sh        Load Geant4 + ROOT from CVMFS
    ├── build.sh               CMake configure + make
    ├── submit_condor_full.py  HTCondor job submission
    ├── plot_geometry.py       2D + 3D detector geometry plots (pyvista)
    └── sample_x17_hits.py     Bridge: Geant4 hit tree → Python fast-MC format
```

---

## Notes

1. **He-3 isotope**: defined with `G4Isotope` (A=3, Z=2), not natural helium.
   Geant4 HP physics includes the ³He(n,p)T thermal capture cross section.

2. **Pair kinematics**: the generator reproduces the Python script
   `e+e-_rel_angle_sim.py` — isotropic decay in the X17 rest frame, then a
   Lorentz boost along +Y with β = p/E at `x17Energy_MeV`.

3. **Thread safety**: each worker thread fills its own `EventData` vector and
   writes its own ROOT file. No locks needed. Merge with `hadd` after all
   threads complete.

4. **LS optical photons**: LAB is defined for dE/dx transport only. Scintillation
   photon propagation requires adding optical material properties and registering
   `G4OpticalPhysics` — not currently enabled.

5. **Step size in He-3**: 1 mm default. Reduce in `SimConfig` if spatial track
   detail inside the target is needed.

6. **Geometry visualisation**: requires the Python venv at
   `~/PycharmProjects/nTof_x17/venv` (has pyvista + matplotlib).
   ```bash
   ~/PycharmProjects/nTof_x17/venv/bin/python scripts/plot_geometry.py
   ```
