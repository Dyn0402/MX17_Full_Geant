# Plan: Geant4 background campaign + Python pile-up layer

Goal: replace the toy IPC-background estimates (`MX17_Simulation/run_ipc_background_study.py`)
with full-simulation numbers, produce all event classes with the new STEP-derived
target geometry, and hand slimmed event pools to the Python fast-MC, which adds
the time structure (pile-up, combinatorics, dead time).

The factorization principle: **Geant4 owns per-interaction physics**
(transport, scattering, capture, conversion, detector response);
**Python owns the time domain** (rates, TOF, coincidence windows, DAQ).
Geant4 events are timeless "interaction templates"; Python composes them
into pulses.

---

## Stage 0 — preparation (local, ~1 day)

1. **Geometry sanity**: `plot_geometry.py` images of the STEP vessel;
   confirm on-axis material budget. Note: the new vessel has
   **~5 mm Al (dome tip) on the beam axis upstream and ~21 mm Al
   (neck + valve) downstream** — ~10–40× more Al in the beam than the old
   flat 0.5 mm end caps assumed in the toy study. Capture rates will rise
   accordingly; this is the single biggest motivation for the neutron run.
2. **Physics validation mini-runs** (1e6 neutrons, local):
   - ³He(n,p) capture rate in gas vs analytic 1/v expectation (G4NDL/HP).
   - Capture-γ spectrum from Al: compare to IAEA-PGAA lines (7.72 MeV etc.).
   - Verify thermal treatment (`G4NeutronHPThermalScattering` not needed for
     gas, but check elastic scattering in Al/CFRP doesn't soften the
     spectrum unphysically).
3. **Known Geant4 limitation — document in README**: Geant4 does **not**
   simulate internal pair conversion (IPC). Photon evaporation emits real γs
   (+ conversion electrons), never internal e⁺e⁻ pairs. Therefore:
   - the neutron run covers **external conversion** of capture γs (which the
     toy says dominates wall pairs ~10:1 over wall IPC) and **all singles**;
   - He-3 IPC and X17 remain **explicit generator modes** (as now);
   - wall IPC (the residual ~10%) can be added later as a generator mode
     reading capture vertices (Stage 2c, optional).

## Stage 1 — new generator modes (C++, ~1–2 days)

Extend `X17PrimaryGenerator` (mode pattern already in place):

1. `--neutron <Emin_eV> <Emax_eV>`: fire neutrons along +Y from upstream of
   the vessel.
   - Energy: iso-lethargic (1/E) between Emin and Emax (default 0.025 eV –
     1 keV); optionally `--nflux <file>` to sample a real n_TOF EAR2 flux
     histogram later.
   - Beam profile: `--beam-fwhm <cm>` Gaussian (or flat disk) — **ask the
     n_TOF beam group for the EAR2 profile**; if the beam is wider than the
     20 mm capsule bore, the halo hitting the thick Al shoulder/neck
     matters a lot.
   - Record per event in EventTree: E_n, vertex of first capture
     (`nCapture` position + volume name), so TOF can be assigned offline.
2. `--gamma-source <volume>`: re-emit capture γs from a capture-vertex
   library (text/ROOT file produced by a neutron run) with the Al/C/H
   cascade energies. This is the **biasing stage** for double-trigger
   statistics — see Stage 3.
3. Keep `-n/--ipc/--single` modes unchanged (X17 + He3-IPC + calibration).

Output schema: unchanged HitTree (+ `event_type` = 2 neutron, 3 gamma-source)
so `analyze_pairs.py` machinery keeps working.

## Stage 2 — production runs (HTCondor, parallel with Stage 3 development)

| Run | Generator | Events | Purpose |
|-----|-----------|--------|---------|
| A | pairs, X17+IPC 50/50 | 10M (100×100k) | signal/IPC pools, response JSON, MS budget for NEW geometry |
| B | `--neutron` | 1e9 (100×1e7) | capture budget, singles rates & spectra, external-conversion background, direct (unbiased) wall double-trigger rate |
| C | `--gamma-source He3Cap_Al` (+CFRP) | 1e8 γ | high-stats wall double-trigger sample (biased; weight = captures/γ from run B) |
| D | `--single e- KE scan` | 10×100k | clean per-particle trigger curves P(trig │ KE, θ) for the response JSON |

Statistics check for run B: toy says ~1.4e-7 wall double-trig/neutron
(old thin caps) → with thicker Al expect ~1e-6–1e-5: 1e9 neutrons give
1e3–1e4 direct events — enough to *validate* run C's biased estimate.
Singles need only ~1e6 captures — trivially covered.

## Stage 3 — analysis + slimming (Python, ~3–4 days)

1. Extend `analyze_pairs.py` (or a sibling `analyze_neutrons.py`) with:
   capture maps (volume × E_n), singles rate per detector per capture,
   conversion-pair spectra, trigger classification — same accumulator
   pattern.
2. **`scripts/make_event_pools.py`** — the slimming step:
   - input: HitTree files from runs A/B/C;
   - output: one compact file per source class
     (`pool_x17.npz`, `pool_ipc_he3.npz`, `pool_wall_bg.npz`,
     `pool_singles.npz`) containing per event: arm/detector IDs, hit
     positions (u,v), times relative to event start, edep digests, truth
     labels (event_type, true KEs, true opening angle, capture E_n);
   - keep it detector-summary-level (one row per fired detector per event),
     ~100 bytes/event, so 10M events ≈ 1 GB → trivially loadable.

## Stage 4 — Python pile-up layer (~3–4 days)

New sampler in `MX17_Simulation` (replaces the response-JSON smearing for
pairs — pool events carry exact correlations):

1. Per pulse: draw N_src ~ Poisson(rate_src) for each pool; assign each
   event a TOF time (from its E_n for neutron-induced classes; from the
   X17/IPC TOF spectrum for signal); inject all hits into the pulse
   timeline.
2. Reuse the existing machinery: `merge_hits`, coincidence pairing,
   trigger logic, plus `dead_time_sim`'s veto/dead-time model — now acting
   on physically correct mixtures.
3. Outputs: combinatorial spectra with real wall-background energies
   (replacing the toy's analytic estimate), trigger rates, dead-time
   losses, and the final profiled significance with every background in.
4. Decision deliverables: Config A vs B re-check, per-leg energy threshold
   optimization (the toy says ~5–6 MeV — verify), beam-window/material
   recommendations.

## Rates normalization

Tie everything to captures: run B gives captures/neutron per volume;
Alberto's chain gives He3-IPC/pulse = 1.12e-2/0.3 — combined with the
He3-IPC branching (54 µb / 5333 b × α_IPC) this fixes neutrons-on-target
per pulse, which then normalizes ALL pools (wall γ, singles, X17 at the
assumed 2.5%). Single consistent ladder, no independent flux guess needed.

## Open questions to settle early

1. **EAR2 beam profile and halo** at the target position (ask beam group) —
   determines how much flux hits the thick Al shoulder directly.
2. Does the beamline have collimation/shielding entering the world volume
   that should be modeled (room-return neutrons are NOT in the sim)?
3. Is the valve (downstream, +Y) really in-beam, or is the capsule mounted
   valve-up perpendicular...? (Vessel axis = beam axis per the STEP
   placement — confirm orientation.)
4. G4NDL completeness for ³He(n,γ): likely absent/negligible — He3-IPC is
   generated explicitly anyway, but confirm no double counting.
