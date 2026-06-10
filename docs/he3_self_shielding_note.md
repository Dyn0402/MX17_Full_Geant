# He-3 self-shielding and the IPC/X17 rate normalization

**Date:** 2026-06-10 · **Status:** needs discussion with Alberto + Geant4 verification

## Summary

The IPC (and therefore X17) rates per pulse in `calculation_tables/results_3He`
appear to be computed with a **thin-target** radiative-capture probability,
`P_rad = (areal density) × σ_nγ(E)`. Below ~1 keV the He-3 gas is **optically
thick** to ³He(n,p)t — the (n,p) optical depth at the tabulated areal density
(3.158×10⁻² at/b) is ~150 at 25 meV — so a sub-keV neutron is absorbed within
a fraction of a millimetre and never samples the full gas column. In that
regime the radiative capture probability per *incident* neutron saturates at
the branching ratio

P_rad(max) = σ_nγ / σ_np = 54 µb / 5333 b = **1.01×10⁻⁸**,

independent of gas thickness. The thin-target formula exceeds this cap by
exactly the optical depth.

## Decomposition (Alberto's table rows vs physical cap)

| E range | nLσ(n,p) | table rad. capt./n | physical cap | overestimate |
|---|---|---|---|---|
| 0.01–0.1 eV | 151 | 1.41×10⁻⁶ | 1.01×10⁻⁸ | ×139 |
| 0.1–1 eV    | 48  | 5.71×10⁻⁷ | 1.01×10⁻⁸ | ×56 |
| 1–10 eV     | 15  | 1.68×10⁻⁷ | 1.01×10⁻⁸ | ×17 |
| 10–100 eV   | 4.8 | 6.69×10⁻⁸ | 1.00×10⁻⁸ | ×6.7 |
| 100–1000 eV | 1.5 | 5.90×10⁻⁸ | 7.9×10⁻⁹  | ×7.5 |

Each table row is reproduced (to ~10%) by
`3.158e-2 at/b × 54 µb × √(25.3 meV/E)` — confirming the thin-target origin.

Flux-weighted over 0.01 eV–1 keV:

- Table IPC rate: **1.21×10⁻² /pulse** (this is the `IPC_PER_PULSE ≈ 1.12e-2`
  used in the sensitivity studies, which corresponds to the E < 1 eV rows)
- Self-shielding-corrected (same α_IPC = 2.1×10⁻³): **1.9×10⁻⁴ /pulse**
- **Net factor: ×62**

Secondary difference: the IPC coefficient (table: 2.1×10⁻³/capture; Rose-table
E1 estimate at 20.6 MeV: ~3.6×10⁻³). Factor ~1.7, choice of multipolarity —
the 2.1×10⁻³ may well be the appropriate value; immaterial next to the ×62.

## Consequences if confirmed

- X17 and He3-IPC scale together by ~1/62 in the sub-keV ROI →
  discovery significance scales by √ → the 3.0σ/30 d projection becomes
  **~0.4σ** if the analysis is confined to E_n < 1 keV.
- **The high-energy region is not affected.** Above ~100 keV the gas is thin
  (optical depth < 0.05) and the thin-target numbers are valid: the table's
  0.1–10 MeV rows carry 6.9×10⁻² IPC/pulse — ~98% of the physical total.
  The commented-out 0.2–2 MeV row at the bottom of `results_3He`
  (5.2×10⁻² IPC/pulse, Δt ≈ 2 µs TOF window) suggests this region was
  already under consideration.
- Moving the ROI to ~0.1–10 MeV would preserve the statistics but changes
  the generator physics: capture at E_n adds energy and a CM boost to the
  4He* system (E* ≈ 20.58 MeV + ¾E_n), shifting the e⁺e⁻ kinematics and
  opening angle. The current at-rest 20.58 MeV generator would need
  updating, and the gamma-flash / TOF veto situation at ~µs timescales
  needs assessment.

## Verification path

The `--neutron` mode of `mx17_full_sim` (this repo) does full neutron
transport with G4NDL/HP data — self-shielding, scattering, and resonances
are automatic. Comparing captures/neutron per energy decade from the
neutron run against both the table and the analytic cap settles this
independently of any analytic assumption. Runs submitted 2026-06-10:
sub-keV window and full-range (1 meV–100 MeV).

Open questions for Alberto:
1. Was the sensitivity anchored to the sub-eV rows (1.12×10⁻²/pulse) or to
   the MeV region?
2. Is the IPC/capture = 2.1×10⁻³ coefficient from an M1 evaluation (and
   should E0 strength be included for this transition)?
3. The table geometry is a 4 cm sphere; the actual vessel is the STEP
   capsule (2 cm bore) — areal densities differ by ~×1.5, also worth
   refreshing.
