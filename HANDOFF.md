# Handoff

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
