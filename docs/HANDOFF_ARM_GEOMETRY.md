# Handoff — investigate the 4-arm geometry vs the as-built module

**Written 2026-08-06**, at the end of the single-module as-built geometry
upgrade. This is a work order for whoever (human or model) picks up the 4-arm
question next. The single-module model is done and validated; do **not**
change `../MX17_Geant/shared/MX17ModuleGeometry.hh` to make the arms fit —
the module is the part backed by CAD + gerbers. The thing to investigate is
the **arm placement description in this repo**.

## The observation

When the as-built module (CAD-derived, 440 mm frame, 470 mm plates offset
(+15,+15) from the active axis, M1 front-end cards straddling the board edges
out to 260.4 mm) is placed at this repo's surveyed arm positions
(`SimConfig.hh`: `mm_distance_x_cm = 20.40`, `mm_distance_z_cm = 20.45`,
pinwheel shifts 1.55 / 1.575 / 1.635 / 1.73 cm), the G4 overlap checker finds
**mm-level interpenetrations between adjacent arms**:

1. **M1 front-end cards ↔ neighbouring arm's window flange**: the cards
   reach 260.4 mm from their arm's active axis and poke ~1 mm into the
   neighbour's flange-ring corner region. (First seen with a wrong M1 model
   at 230.8 mm reach — the corrected, gerber-aligned cards reach further, so
   the clash is worse, ~10 mm scale.)
2. **Window bulge ↔ neighbouring flange corner**: the terraced-dome model is
   square, so its outer terraces carry the full aperture footprint at 2 mm
   height; the corners would graze the neighbour's flange. This one is at
   least partly a model artifact (a real dome falls to zero at the clamped
   edge), but the margin is sub-mm either way.

**Dylan reports the real assembly in the experimental area has no
interpenetration.** Since the module dimensions are solid (CAD + two
independent gerber cross-checks), *some part of the arm-placement description
here is probably wrong or incomplete.*

Current state of this repo pending that investigation (both commented in
`src/DetectorConstruction.cc`): `mmSpec.feThick_mm = 0` (M1 cards off) and
`AsBuiltSpec(0.0)` (bulge off).

## Hypotheses, most likely first

1. **Per-arm orientation of the connector edges is wrong.** The shared module
   places the M1 cards and the (+15,+15) plate offset on the module-local
   +x/+y edges, and every arm currently uses the same local orientation. If
   the real modules are mounted so the connector edges face *outward* (away
   from the neighbouring arm) — which is exactly what a designer would do —
   the cards land in free space and the clash disappears. This is a pure
   bookkeeping fix: a per-arm rotation (0/90/180/270°) of the module about
   its normal. **Check first**: photos of the assembled station, or the cable
   routing (M1 cards are where the VMM cables attach — cables visibly exit on
   specific sides of each arm).
2. **Arm front-face distances.** 204.0 / 204.5 mm were measured to the
   drift-mylar front face (2026-07-15, pre-as-built). If they were actually
   measured to the *window flange* or *frame front*, the true window-plane
   distance is up to ~5 mm larger, which relieves the corners.
3. **Pinwheel shifts.** 15.5–17.3 mm measured per arm; the clash regions are
   at the arm corners, where a ±5 mm shift error changes clearances one-for-one.
4. **The flange/frame outer profile is not a full 440 mm square.** The
   gas-frame gerber outline (`MX17_Geant/design/gerbers/readout_pcb/
   DFS3498A_gasframe.gbr`) has an asymmetric outline (bbox 451 × 447,
   centred ~(+7,+9)) with tabs — the sim's plain 440 ring may be too fat
   exactly at the corners where the neighbour's M1 cards pass.

## Concrete tasks

1. Determine the per-arm module orientation (hypothesis 1) from photos /
   cable exits; add a per-arm quarter-turn parameter to the module placement
   in `src/DetectorConstruction.cc` (the shared header needs no change — the
   pieces are placed by the consumer, so a rotation about the module normal
   in the placement loop is ~10 lines).
2. Re-enable `feThick_mm` (M1 cards) and re-run the overlap check per
   orientation hypothesis; the checker is definitive and cheap
   (`mx17_full_sim -n 1 -t 1 --single electron 5 90 0`).
3. If clashes persist, re-measure/confirm the arm distances and pinwheel
   shifts against hypothesis 2/3, or model the true frame outline
   (hypothesis 4) from the gasframe gerber.
4. Re-enable the window bulge once clearances are understood. If the corner
   grazing is the only remaining issue, either accept the flat window in this
   repo (documented) or implement a rounded-corner terrace profile
   (`ScaleAbout`-style shrink of a rounded-square outline) in the shared
   header.

## Reference numbers

Module footprints, all referenced to the module's ACTIVE-AREA axis
(see `MX17_Geant/design/figures/mx17_plan_views.png`):

| element | extent |
|---|---|
| window / flange / frame outer | 440 mm square, centred |
| 470 mm plates (board, rohacell, support plate) | −220..+250 in x and y |
| M1 cards (2 on +x edge, 2 on +y edge) | radial 219.4..260.4, tangential ±100 ± 80 |
| window bulge footprint | 400.2 mm square, sag 8 mm upstream |
| arm front faces (this repo, survey) | 204.0 mm (±X), 204.5 mm (±Z) |
| pinwheel shifts (arm 0..3) | 15.5 / 15.75 / 16.35 / 17.3 mm along −û |

Validation harness used during the module work (fixed-seed digests,
overlap-check greps): see `MX17_Geant` commit 54f3544 and the notes in
`MX17_Geant/design/GEOMETRY_IMPLEMENTATION_NOTES.md`.
