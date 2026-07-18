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

## 5. Detector stack, inside → out (arrangement flip, 2026-07-15)

The layer order was flipped and the liquid reduced to a single layer.  New
order: **MM → SiPM wall → plastics → 1 LS layer.**  All depths are quoted from
the **MM drift-mylar front face** (w = 0).  The MM/PCB build is unchanged; the
air gap after the PCB simply grows so the SiPM container lands at the measured
distance.

| Layer | Depth (front → back) | Centered on | Notes |
|-------|----------------------|-------------|-------|
| MM + PCB | 0 → 3.60 cm | MM (pinwheel) | unchanged internal build |
| **SiPM wall** | **11.0 → 14.5 cm** (container 3.5 cm, measured 2026-07-17; scint centered) | **STRUCTURE** (u=0) | 50×50 cm, 20 bars × 2.5 cm; **16 read out**, window shifted 1 bar toward the MM (drop 3 far, 1 near); un-read bars removed in sim / transparent in plots |
| Plastics | per arm (see below) | MM (pinwheel) | 2×(20×30 cm) bars, **2.5 cm** thick (nominal 2.0) |
| LS vessel | per arm (see below) | MM (pinwheel) in u — **u NOT measured, assumed** | STEP-derived shape; per-arm depth/orientation/height surveyed 2026-07-17 |

### Per-arm survey 2026-07-17 (all from the SiPM container back = 14.5 cm)

| Arm | Wall | SiPM back → plastics front | Plastics (abs) | SiPM back → LS flat slab face | LS slab (abs; apex ±1.25 cm more) | LS orientation | LS slab-centre v |
|-----|------|---------------------------|----------------|-------------------------------|------------------------------------|----------------|------------------|
| 0 | D (+X) | 6.5 cm | 21.0 → 23.5 cm | 12.3 cm | 26.8 → 28.9 cm | HORIZONTAL, PMT +u (−Z, West) | −0.04 cm |
| 1 | B (−X) | 6.1 cm | 20.6 → 23.1 cm | 12.7 cm | 27.2 → 29.3 cm | VERTICAL, PMT up (+v) | +0.03 cm |
| 2 | A (+Z) | 6.3 cm | 20.8 → 23.3 cm | 12.3 cm | 26.8 → 28.9 cm | HORIZONTAL, PMT +u (+X, North) | +0.06 cm |
| 3 | C (−Z) | 6.1 cm | 20.6 → 23.1 cm | 12.7 cm | 27.2 → 29.3 cm | VERTICAL, PMT up (+v) | −0.07 cm |

LS depth reference = the **flat slab front face** (measured at the vessel edge,
off the bulge); the front bulge apex sits 1.25 cm (=hCap) closer to the target.
"PMT +u" = "to the right when looking from behind the wall toward the target"
(up = sky).  Heights from the bottom-bar chain: SiPM enclosures 62 cm tall
(SiPM bars assumed centred → enclosure bottom at v = −31); a bar sits
**6.8 cm** above the enclosure bottom (**DOUBLE-CHECK this 6.8**); LS vessel
bottoms above that bar: D 1.6, B 1.7, A 1.7, C 1.6 cm.  With slab half-heights
22.53 (vertical) / 22.56 cm (horizontal), all four slab centres land within
±0.7 mm of v = 0 — a strong consistency check.  B/C PMT-up is *derived* from
this chain (PMT-down would put the slabs ~27 cm off-centre).

### LS vessel (STEP-derived, 2026-07-17)

Source: `LS X17.step` (Shapr3D export 2026-06-10, `/media/dylan/data/x17/LS_Stuff/`;
single solid 'Corpo 03', products SCINT1/SCINT2).  Modelled in
`src/DetectorConstruction.cc` + `scripts/plot_geometry.py` (kept in sync):

- **Slab**: 45.12 × 45.06 cm face, 21.2 mm outer thickness; CFRP wall 2.6 mm
  (2.0 structural + 0.6 liner; the 40 µm Al liner is neglected).
- **Funnel**: 9 cm loft from the slab cross-section down to the Ø50 mm neck
  (STEP B-spline loft approximated as a G4Trd; square→circle mismatch at the
  neck end is negligible).
- **Neck**: r = 25 mm × 12.97 cm; **PMT inserted halfway** (window 6.5 cm into
  the neck, r = 2.2 cm borosilicate envelope + vacuum, ~11.5 cm long).
- **Liquid (scored `LiqScint_1`)**: LAB filling slab + funnel + neck up to the
  PMT window.  The **6.5 L fill** exceeds the flat interior (~3.9 L), so the
  slab faces **bulge** — modelled as two ellipsoidal domes whose height is
  solved from the fill volume: **12.5 mm per side** (volume-exact; footprint
  is elliptical instead of the true rounded-square pillow).
- **Orientation** (surveyed 2026-07-17): VERTICAL with PMT up on B/C,
  HORIZONTAL with PMT along +u on A/D (`ls_rot_deg` in `SimConfig.hh`).
- Dropped STEP details: r=3.5 mm fill tube along the face, r=9 mm boss,
  r=4–10 mm edge fillets.

Measured values: SiPM container front 11 cm from mylar front, container depth
3.3 cm.  Gaps flagged "measure later" are rough.  See
`GEOMETRY_CHANGE_CHECKLIST.md` for the full list of files that must move together.

## 6. Propagation status

**Propagated into the code 2026-07-14 (coord convention + pinwheel):**

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

**Propagated into the code 2026-07-15 (stack flip; see §5):**

- [x] `include/SimConfig.hh` — new SiPM/plastics/LS-gap fields; removed the old
      `scint_size_*` and `gap_*` fields; plastics thickness 2.0→2.5 cm.
- [x] `src/DetectorConstruction.cc` — SiPM wall = 16 bars on the structure,
      plastics + single LS box on the MM; per-axis absolute depths.
- [x] `include/DetectorConstruction.hh` / `src/SteppingAction.cc` — dropped the
      2nd LS layer (`fLS2LV` / `LiqScint_2`).
- [x] `scripts/plot_geometry.py`, `plot_buildup.py` — new stack, structure-vs-MM
      centering, transparent un-read SiPM bars; 2D top-down now uses the adopted
      +Z-right/+X-up orientation (legacy-orientation note resolved).

**Propagated into the code 2026-07-17 (STEP LS vessel; see §5):**

- [x] `include/SimConfig.hh` — box-LS fields replaced by STEP vessel fields
      (`ls_slab_*`, `ls_wall_mm`, `ls_funnel_len_cm`, `ls_neck_*`,
      `ls_fill_liters`, `ls_pmt_*`, `ls_neck_dir`).
- [x] `src/DetectorConstruction.cc` — CFRP shell union (slab+funnel+neck+2
      bulge domes), LAB liquid daughter, PMT bore/glass/vacuum; overlap checks
      clean; hit u/v ranges verified (liquid reaches the PMT window at
      v = +380 mm, +v only).
- [x] `src/SteppingAction.cc` — `LiqScint_1` arm ID now walked 1 level up
      (nested in the vessel shell).
- [x] `scripts/plot_geometry.py` / `plot_buildup.py` — vessel outlines (bulge
      lens in top-down, full silhouette in side view, 3D vessel); regenerated.

**Propagated into the code 2026-07-17 evening (placement survey; see §5):**

- [x] `include/SimConfig.hh` — SiPM container 3.3→3.5 cm; per-arm
      `gap_sipm_to_plastic_cm[4]`, `ls_front_from_sipm_back_cm[4]`,
      `ls_rot_deg[4]`, `ls_offset_v_cm[4]`; removed `gap_plastic_to_ls_cm`
      and `ls_neck_dir`.
- [x] `src/DetectorConstruction.cc` — per-arm plastics/LS depths, vessel
      rotation about w (vertical B/C, horizontal A/D), surveyed v offsets;
      overlap checks clean; per-arm hit u/v ranges verified.
- [x] `scripts/plot_geometry.py` / `plot_buildup.py` — per-arm depths,
      orientation-aware LS drawing; all figures regenerated.

**Still open:**

- [ ] Re-run simulations & acceptance with the new geometry — existing sim
      outputs are now stale (again after the 2026-07-17 changes).
- [ ] **Double-check the 6.8 cm bar height** used as the LS-bottom reference.
- [ ] Measure the LS horizontal positions along each wall (u) — slab
      horizontally-centred-on-MM is assumed; for the horizontal vessels (A/D)
      this also fixes where the funnel/PMT sit.  (E.g. measure SiPM-enclosure
      edge → LS slab edge per wall.)  NB: the bottom-bar survey fixed the
      VERTICAL (v) positions; u is the left–right direction along the wall.
- [x] PMT directions CONFIRMED 2026-07-18: D→West/−Z, A→North/+X, B/C up.
- [ ] `sipm_front_from_mylar_cm = 11.0` not re-measured in this survey.
- [ ] Trim analysis scripts that reference `LiqScint_2` (harmless: no hits now).
- [ ] Measure & update the "measure later" gaps (SiPM→plastics 7 cm,
      plastics→LS 5 cm) and the LS MM-centering.
- [ ] Slides / reports quoting the old stack (25 cm front face, 2 LS layers).

Note: values were 41 cm / uniform 30 mm in the first pass (2026-06-30 early),
then 3.10–3.46 cm shifts; the current numbers above supersede them.
