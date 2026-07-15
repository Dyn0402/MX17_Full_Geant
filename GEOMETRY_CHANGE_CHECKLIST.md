# MX17 Geometry-Change Checklist

**Every time the detector geometry changes, touch each item below.** The source
of truth is `include/SimConfig.hh`; everything else must be brought back in sync
with it. (See `GEOMETRY_COORDINATE_CONVENTION.md` for the coordinate convention
and the current measured numbers.)

| # | File | What lives here |
|---|------|-----------------|
| 1 | `include/SimConfig.hh` | **Source of truth** — all distances, layer sizes, thicknesses, gaps, pinwheel shifts. |
| 2 | `src/DetectorConstruction.cc` | Layer thicknesses, slab lists, per-arm placement (MM / SiPM / plastics / LS), sensitive-detector registration. |
| 3 | `include/DetectorConstruction.hh` | Sensitive-volume member pointers (add/remove when a scored layer is added/removed). |
| 4 | `src/SteppingAction.cc` | `kScoredVolumes` — the set of scored volume *names*. Update if a scored volume is renamed/added/removed. |
| 5 | `mx17_full_sim.cc` | CLI overrides that write geometry fields (e.g. `--dist`). |
| 6 | `scripts/plot_geometry.py` | `CFG` dict (mirror of SimConfig) + the layer-depth chain + `LAYERS_2D` + the 3D `layer_info`. Top-down / side / 3D figures. |
| 7 | `scripts/plot_buildup.py` | Imports depths from `plot_geometry`, but has its OWN `MMS`/`SPAN`, stage names, and `_draw_arm_layers`. Inside-out build-up slides. |
| 8 | `scripts/plot_mm_layout.py` | MM-only canonical top-down. Only needs edits if the **MM** geometry (spacing / pinwheel / size) changes. |
| 9 | `GEOMETRY_COORDINATE_CONVENTION.md` | Human record: measured numbers, tables, propagation status. |
| 10 | Memory: `…/memory/coord-convention-mm-layout.md` + `MEMORY.md` | One-line + detail memory of the adopted convention & current numbers. |
| 11 | `~/PycharmProjects/nTof_x17_DAQ/run_config_beam.py` (separate repo) | Detector *positions* for the DAQ. Only needs edits if the numbers it stores (currently MM mesh-plane positions) change. The 2026-07-15 SiPM/plastics/LS change did **not** move the MM, so DAQ was untouched — re-check if MM positions ever change. |
| 12 | Slides / reports (`.claude/**`, slide decks) | Any figure or number quoting the geometry. Flagged, not auto-updated. |

## Downstream (re-run, not code edits)

- **Sims are stale** after any geometry change — existing outputs used the old
  geometry. Re-run + re-run acceptance/analysis.
- Analysis scripts that hard-code scored volume names
  (`analyze_pairs.py`, `check_output.py`, `make_event_pools.py`,
  `analyze_mev_captures.py`, `analyze_thermal_captures.py`,
  `make_thermal_report_figures.py`, `study_invmass.py`) reference `LiqScint_2`.
  With the single LS layer that volume simply has no hits — harmless, but those
  scripts should be trimmed to one LS layer when analysis is next revisited.

## Regenerate figures after a change

```bash
source ../../PycharmProjects/nTof_x17/.venv/bin/activate   # (see memory)
python scripts/plot_mm_layout.py
python scripts/plot_geometry.py            # add --no-3d to skip pyvista
python scripts/plot_buildup.py             # both styles, 4 frames each
```
