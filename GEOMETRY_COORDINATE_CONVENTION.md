# MX17 Coordinate Convention & MM Layout

**Status:** adopted 2026-06-30 at the experiment. The coordinate *convention*
below matches the existing Geant simulation. The **MM spacing/shift measurements
are new and NOT yet propagated** to the sim / `SimConfig` / analysis — they must
be fixed everywhere later (see "To fix elsewhere" at the bottom).

Canonical drawing: `scripts/plot_mm_layout.py` → `scripts/mx17_mm_layout_topdown.png/.pdf`.
Inside-out build-up slides (MM → +SiPM → +LS → +PMT, two styles):
`scripts/plot_buildup.py` → `scripts/mx17_buildup_{clean,detailed}_{1_mm,2_sipm,3_ls,4_full}.png/.pdf`.

---

## 1. Coordinate system (right-handed)

- **Beam = +Y**, floor → ceiling (`X17PrimaryGenerator.cc`:
  `SetParticleMomentumDirection(0,1,0)`; `DetectorConstruction.hh:4`).
- Right-handed: x̂ × ŷ = ẑ.
- Origin at the He-3 target centre.

### Top-down drawing convention
A **true top-down** view (looking straight down from above):

| Axis | Screen direction | Cardinal |
|------|------------------|----------|
| +X   | up               | **North** |
| +Z   | right            | **East**  |
| +Y (beam) | out of the page (⊙, toward viewer) | — |

Screen mapping used in the plot: `screen_x = Z`, `screen_y = X`.

Why this works: with +Z→right and +X→up, right-handedness gives
ŷ = ẑ × x̂ = (right) × (up) = **out of the page**, so +Y (beam) points at the
viewer — a genuine bird's-eye view that is still consistent with the sim's +Y
beam. The map (sim X→up, Y→out, Z→right) is a cyclic permutation of axes = a
**pure rotation** (determinant +1, no mirror), so sim and diagram stay
compatible.

> Note: you cannot simultaneously have +X-right, +Z-up, beam-out-of-page in a
> *left*-handed reading. The earlier "+X right / +Z up" attempt forced beam=−Y;
> swapping to **+Z right / +X up** restores beam=+Y. This swap is the adopted
> convention.

---

## 2. Detector labels (A / B / C / D)

| MM | Coord | Screen position | Cardinal | Sim arm |
|----|-------|-----------------|----------|---------|
| **A** | +Z | right  | East  | Arm 2 (+Z) |
| **B** | −X | bottom | South | Arm 1 (−X) |
| **C** | −Z | left   | West  | Arm 3 (−Z) |
| **D** | +X | top    | North | Arm 0 (+X) |

(Sim arm IDs from `DetectorConstruction.cc:515-520`.)

---

## 3. MM measurements (refined 2026-06-30)

Opposing **mylar (entrance-window) face** distances:

| Pair | Distance | Faces at (centred) |
|------|----------|--------------------|
| **B (−X) ↔ D (+X)** | **40.8 cm** | ±20.40 cm in X |
| **C (−Z) ↔ A (+Z)** | **40.9 cm** | ±20.45 cm in Z |

| Quantity | Value | Notes |
|----------|-------|-------|
| MM active width (in-plane, u) | 38 cm | unchanged |
| MM stack depth | ≈ 3.04 cm | incl. 30 mm drift gap |

**Recentering:** the set-up is positioned so the **target/beam sit at the centre
of the (roughly square) box formed by the four inner mylar faces**. Because the
square's centre is the midpoint of each opposing face pair, this means each
mylar face is half its pair-span from the beam axis — X faces at ±20.40 cm, Z
faces at ±20.45 cm. (The asymmetric tangential shifts do not move the faces
perpendicular to themselves, so they don't affect centering.)

---

## 4. Pinwheel (circular) shift — per detector

Every MM is shifted **tangentially** (⟂ its outward normal) in a right-handed
pinwheel. For a detector with outward normal **n̂ = (n_x, n_z)** (physics X,Z):

```
shift = t̂ · |shift|,   where   t̂ = (−n_z, n_x)   (in-plane tangent)
```

Magnitudes are **individual** (no longer a uniform 30 mm):

| MM | Position | Shift direction | **Magnitude** |
|----|----------|-----------------|---------------|
| **D** | +X (top)    | → right (+Z) | **1.55 cm** |
| **A** | +Z (right)  | → down (−X)  | **1.635 cm** |
| **B** | −X (bottom) | → left (−Z)  | **1.575 cm** |
| **C** | −Z (left)   | → up (+X)    | **1.73 cm** |

i.e. a clockwise pinwheel in the top-down view. (User-specified examples:
"B shifted left, C shifted up" — both reproduced by the formula.)

> **Corrected 2026-07-14:** the earlier magnitudes (D=3.10, A=3.27, B=3.15,
> C=3.46 cm) were **2× too large**; the true measured offsets are half those
> values, tabulated above.

---

## 5. Propagation status

**Propagated into the code 2026-07-14:**

- [x] `include/SimConfig.hh` — replaced `mm_distance_cm=25` with per-axis
      `mm_distance_x_cm=20.40` / `mm_distance_z_cm=20.45`, added
      `mm_pinwheel_shift_cm[4]={1.55,1.575,1.635,1.73}` (arm order D,B,A,C).
- [x] `src/DetectorConstruction.cc` — per-axis front-face distances in `armDefs`;
      each arm slid along −uHat by its pinwheel amount (the whole assembly,
      incl. its local coordinate origin → hit u-coords stay centred). Verified:
      arm front faces print (20.4,0,1.55) / (−20.4,0,−1.575) / (−1.635,0,20.45) /
      (1.73,0,−20.45) cm.
- [x] `mx17_full_sim.cc` — `--dist` now sets both per-axis distances (uniform override).
- [x] `scripts/plot_geometry.py` — `CFG` per-axis + pinwheel; `ARM_DEF` builds
      shifted front faces; 2D top-down / side / **3D** figures regenerated.
- [x] `include/RunAction.hh` — added missing `#include <fstream>` (latent bug
      surfaced when the SimConfig change forced a recompile; unrelated to geometry).

**Still open:**

- [ ] Re-run simulations & acceptance with the new geometry — existing sim
      outputs used the old 25 cm symmetric placement and are now stale.
- [ ] Slides / reports quoting the 25 cm front face or symmetric arm placement.
- [ ] `plot_geometry.py`'s 2D top-down still uses its legacy orientation
      (X→right, Z→up, beam into page), not the adopted +Z-right/+X-up view —
      geometry is correct, only the on-page orientation differs from the
      canonical `plot_mm_layout.py` / `plot_buildup.py`.

Note: values were 41 cm / uniform 30 mm in the first pass (2026-06-30 early),
then 3.10–3.46 cm shifts; the current numbers above supersede them.
