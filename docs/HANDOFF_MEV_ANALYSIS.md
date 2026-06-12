# Handoff: MeV-region analysis (start here)

**Written:** 2026-06-12, end of the thermal-statistics session.
**For:** a fresh session analyzing the high-energy (MeV) region from run
B-full, plus the pure X17/IPC pair sample.

---

## Decision context (why this session exists)

The thermal/sub-keV X17 measurement is **dead and the collaboration agrees**
(discussed with Alberto 2026-06-11/12). Reason: below 1 keV the ³He gas is
opaque to (n,p)t, capping radiative captures at σ_nγ/σ_np ≈ 1×10⁻⁸ per beam
neutron — measured with 10⁹ fully-transported Geant4 neutrons as
**0.06–0.11 X17 produced per day** (before acceptance), factor 50–100 below
the collaboration rate table. Full story: [report/thermal_note.pdf](report/thermal_note.pdf)
(12 pp, final) and [thermal_campaign_handoff.md](thermal_campaign_handoff.md).

**The experiment must move to the high-energy region (~0.1–10 MeV)**, where
the gas is thin (opacity < 5% above ~100 keV), the thin-target table rows
are valid, and ~98% of the physical radiative-capture rate lives. The
commented-out 0.2–2 MeV row in `/home/dylan/x17/calculation_tables/results_3He`
(5.2×10⁻² IPC/pulse, Δt ≈ 2 µs TOF window) is the candidate anchor.

## The task now

1. **MeV-region capture budget from run B-full** (`neutrons_fullrange/`,
   5×10⁸ events, 1 meV–100 MeV): radiative captures + (n,p) + wall captures
   per energy decade/bin, per-pulse and per-day normalisation. This is the
   high-energy analogue of the thermal scan — `scripts/analyze_thermal_captures.py`
   is the template (its hardcoded sub-keV binning, LOGE −3..3, needs
   extending to ~+8; the Alberto-table overlays need the >100 keV rows).
   Key physics question: confirm the thin-target rows where they claim
   validity, and map where opacity transitions (expect ~10–100 keV).
2. **Pairs analysis** (`pairs_v2_step_target/`, 10⁷ X17+IPC 50/50, new STEP
   geometry): acceptance, trigger rates, response JSON, MS budget —
   `scripts/analyze_pairs.py` (Q1–Q9) already does this; needs running on
   the new sample and `--export-response geant4_response.json` regenerated.
   **Caveat:** these pairs were generated *at rest* (E* = 20.58 MeV,
   thermal kinematics). For the MeV ROI the ⁴He* gets E* ≈ 20.58 MeV + ¾E_n
   plus a CM boost — kinematics change (opening angle, energy sharing,
   TOF). The at-rest sample is still the right baseline for
   detector-response machinery, but a generator extension + new run A will
   be needed for final MeV-ROI projections.
3. **Gamma-flash / TOF window assessment**: MeV neutrons arrive ~µs after
   the flash at EAR2 (~20 m flight path). E_n ↔ TOF mapping per event is
   available offline via `neutron_E_eV`. Quantify what the 0.2–2 MeV
   window looks like in time and what backgrounds share it.

## Data (all validated: every tree readable, exact entry counts)

`/eos/experiment/ntof/data/x17/full_sim/` (n_TOF experiment EOS):

| dir | files | events | generator |
|---|---|---|---|
| `neutrons_fullrange/` | 100 | 5×10⁸ | `--neutron`, E = 1 meV–100 MeV, EAR2 flux+footprint |
| `neutrons_subkev/` | 100 | 10⁹ | `--neutron`, E = 1 meV–1 keV (analysed, done) |
| `pairs_v2_step_target/` | 100 | 10⁷ | `--ipc 0.5`, X17+IPC 50/50, at-rest kinematics |

EventTree branches for neutron runs: `neutron_E_eV`, `capture_vol`,
`capture_proc`, `cap_x/y/z` + the standard pair-truth branches.
HitTree: per-step hits in sensitive volumes (see README "Output").

## Critical pitfalls (cost us time last session)

- **³He(n,p)t is `neutronInelastic`, NOT `nCapture`** in Geant4 HP.
  `nCapture` in `He3Gas` = radiative capture = the X17 channel.
- **X17 branching: 2.5% is per IPC pair** (X17 = captures × 2.1×10⁻³ ×
  0.025), not per capture — a factor-500 trap (fixed commit 8d3453e).
- **Volume name is `He3Cap_Al`** for the vessel, `He3Gas` for the gas.
- **uproot string branches**: convert with `np.asarray(..., dtype=object)`
  immediately (see `docs/uproot_awkward_concat_pitfall.md`).
- **Always validate output files** (`uproot` + `num_entries`) after any
  campaign — the 06-10 quota incident produced GB-sized unreadable files.
- Normalisation anchors: **7.31×10⁶ n/pulse < 1 keV** (sub-keV runs).
  The fullrange run sampled the same flux file over the full range — the
  n/pulse anchor for that window must be re-read from
  `data/fluxEAR2-Ph3_in_different_units.root`
  (hist `flux_n_pulse_NOisolet_100bpd`); do NOT reuse 7.31×10⁶.
  Pulses/day: **1.929×10⁴** (7×10¹² ppp, from `results_3He`).

## Infrastructure

- Repo: `MX17_Full_Geant` (GitHub Dyn0402, main). lxplus clone:
  `/afs/cern.ch/work/d/dneff/git/MX17_Full_Geant` (built, Geant4 11.2/LCG;
  `source scripts/setup_lxplus.sh` for python3+uproot+pdflatex env).
  SSH: `ssh lxplus` works from this machine (GSSAPI).
- Scan template: `scripts/analyze_thermal_captures.py` (multiprocess,
  skip-bad-files, ~15 min/10⁹ events with 6 workers on lxplus).
- Figures/report toolchain: `scripts/make_thermal_report_figures.py`,
  `scripts/plot_he3_xs.py` (ENDF σ from `data/He3.h5`),
  `scripts/plot_geometry.py` (STEP capsule, in sync with
  `DetectorConstruction.cc`), LaTeX note pattern in `docs/report/`.
- Pairs analysis: `scripts/analyze_pairs.py` (Q1–Q9 + response export),
  `scripts/check_output.py` (single-file sanity).
- Condor submission: `scripts/submit_neutrons.py`, `scripts/submit_pairs.py`
  (defaults already point at ntof EOS).

## Open items carried over

- Alberto follow-ups: IPC coefficient provenance (2.1×10⁻³ —
  multipolarity? E0?); exact X17/IPC = 2.5% basis; capsule geometry update
  of the table (it assumed a 4 cm sphere).
- Generator extension for E_n-dependent kinematics (PLAN item; needed
  before final MeV run A).
- Run C (gamma-source wall backgrounds) and run D (single-particle trigger
  curves) not yet launched — sequencing after the MeV budget makes sense.
- PLAN stages 3–4 (event pools + Python pile-up sampler) untouched.
