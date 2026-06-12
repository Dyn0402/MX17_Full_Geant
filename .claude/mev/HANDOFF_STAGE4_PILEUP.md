# Handoff: Stage 4 pile-up sampler (start here)

**Written:** 2026-06-12, end of the MeV-analysis + pairs session.
**For:** a fresh session building the Python pile-up/timing simulation in
`nTof_x17/MX17_Simulation`, consuming the event pools produced today.

---

## Decision context (why this session exists)

Two campaigns settled the physics this week:

1. **Thermal (sub-keV) is dead** — gas opaque to (n,p)t, 0.06–0.11 X17/day.
   `docs/report/thermal_note.pdf`, `.claude/thermal/`.
2. **MeV region is confirmed feasible** (this session, run B-full, 5×10⁸
   events): **714 direct ³He(n,γ) events → 32.7 X17 produced/day** full
   range, **22.5 ± 1.0/day in the 0.2–2 MeV window**, within 10% of the
   thin-target table row. `docs/report/mev_note.pdf`,
   `analysis/mev/mev_rates.json`.

Pairs acceptance (10⁷ at-rest X17+IPC events, new STEP geometry, also this
session): X17 double-trigger **12.4%**, IPC **3.7%** → **~2.8 double-
triggered X17/day in-window** (~84 per 30-day run), trigger-level S/B ≈ 0.09
before mass/angle cuts. `analysis/pairs_v2/pair_analysis_v2.pdf`.

**The goal now** (user's words): use slimmed event pools small enough to
download locally, pull individual events in Python, and simulate the pulse
with a time component — inject neutron-beam backgrounds + signal into a
common timeline to model the time-compact gamma flash environment.

## What was done today (all on main, local + lxplus in sync)

| deliverable | where |
|---|---|
| MeV capture budget note (8 pp) | `docs/report/mev_note.pdf` (+ `.tex`, figs) |
| Per-decade rates + window integrals | `analysis/mev/mev_rates.json` |
| Raw scan histograms + 714 (n,γ) energies | `analysis/mev/mev_captures.npz` |
| Scan + figure scripts (corrected physics) | `scripts/analyze_mev_captures.py`, `scripts/make_mev_report_figures.py` |
| Pairs analysis PDF (17 sections) | `analysis/pairs_v2/pair_analysis_v2.pdf` |
| Detector-response JSON (10⁷ events) | `analysis/pairs_v2/geant4_response.json` — **installed** in `nTof_x17/MX17_Simulation/geant4_response.json` (verified `Geant4Response` loads it; old kept as `geant4_response_old_geom.json`) |
| Pool builder (PLAN Stage 3) | `scripts/make_event_pools.py` |
| Event pools | `analysis/pools/` on lxplus — **still building, check status first** (see below) |

## IN FLIGHT at handoff: pool builds on lxplus

Launched ~13:30 CET as a sequential background job
(`/afs/cern.ch/work/d/dneff/git/MX17_Full_Geant`, PID 3192012):

1. pairs pools (100 files → `pool_x17.npz` + `pool_ipc.npz`, 500k
   events/class subsample) — 24/100 files at 14:05, ETA ~16:00 CET
2. then neutron-bg pool (`pool_neutron_bg.npz` from `neutrons_fullrange/`,
   `--t-max-us 50`)
3. writes `analysis/pools/pools.status` (line `ALL_POOLS_DONE`) when both done

**First action of next session:**
```bash
ssh lxplus "cat /afs/cern.ch/work/d/dneff/git/MX17_Full_Geant/analysis/pools/pools.status 2>/dev/null; \
            ls -lh /afs/cern.ch/work/d/dneff/git/MX17_Full_Geant/analysis/pools/; \
            tail -3 /afs/cern.ch/work/d/dneff/git/MX17_Full_Geant/analysis/pools/*.log"
# then pull them local:
scp "lxplus:/afs/cern.ch/work/d/dneff/git/MX17_Full_Geant/analysis/pools/pool_*.npz" analysis/pools/
```
If it died (no status file, logs stalled): rerun the two
`make_event_pools.py` commands from the docstring of that script — they are
idempotent. Expected sizes: ~165 MB per pairs pool (from 13 MB/39k-event
smoke test), neutron pool unknown (first measurement — could need the
`--t-max-us` window tightened or loosened; check `meta["rate_per_pulse"]`).

## Pool format (what the sampler consumes)

`np.savez_compressed`; `meta` key = JSON string. Per event: truth arrays
(`evt_*`); per hit-digest: flat arrays (`dig_*`) indexed CSR-style —
event i's digests are rows `dig_start[i]:dig_start[i+1]`.

- Digest = one row per (event, trackClass, arm, detector): `edep_MeV`
  (summed), `t_first_ns`/`t_last_ns` (Geant4 global, t=0 at primary launch),
  `ke_first_MeV`, entry/exit (u,v,w). trackClass: 0 = primary e⁻,
  1 = primary e⁺, 2 = secondaries (summed per volume).
- Pairs pools: truth = vertex, em/ep KE + direction unit vectors,
  openingAngle, inv_mass. Generated at t=0, **no E_n** (at-rest sample).
- Neutron pool: truth = `neutron_E_eV`, capture vol/proc (integer-coded,
  codes in meta), `cap_y`, and **`t0_pulse_ns` = TOF(E_n) over 19.3 m**
  (19.5 m flight path − 20 cm gun offset; gun launches at y=−20 cm, t=0).
  Hit pulse time = `t0_pulse_ns + dig_t_*_ns`. Normalisation in meta:
  `rate_per_pulse = n_pool × n_per_pulse / n_events_simulated`,
  n_per_pulse = 2.2628×10⁷ (flux-file integral — NOT 7.31e6, NOT 3.29e7).

## The task: Stage 4 sampler in `nTof_x17/MX17_Simulation`

Per pulse:
1. **Backgrounds:** draw N ~ Poisson(`rate_per_pulse`) from
   `pool_neutron_bg`, sample events with replacement; each carries its own
   `t0_pulse_ns` — the TOF spectrum and the Geant4 thermalization delays
   come for free.
2. **Signal/IPC:** draw event count from the capture rate (window row of
   `analysis/mev/mev_rates.json`: 22.2 capt/pulse in 0.2–2 MeV ×
   α_IPC = 2.1e-3 for pairs); assign each pair an E_n sampled from the
   per-decade capture curve (or the 714 raw energies in `mev_captures.npz`)
   → TOF → injection time; attach a pool event's hits at that time.
3. **Timeline:** merge all digests into one pulse stream; reuse existing
   `MX17_Simulation` machinery (merge_hits, coincidence pairing, trigger
   logic, `dead_time_sim` veto model).
4. **Gamma flash:** NOT in any Geant4 run. Model as t=0 marker + detector
   blind/recovery-time parameter (scan it — the answer must come from beam
   data). Window of interest sits 1.0–3.2 µs after flash.
5. **Deliverables:** trigger rates and dead-time losses vs flash-recovery
   assumption, combinatorial spectra with real wall-background energies,
   in-window pile-up probability, refreshed significance.

In-window per-pulse context (from `mev_rates.json` window block):
3.0×10⁵ (n,p)t in gas (no detector hits — in-gas heat), 8.8×10⁴ LS
H-captures (2.2 MeV γ inside the calorimeter — the dominant detector
background), 2.6×10³ plastic-veto captures, ~255 Al-wall captures vs
22 IPC-capable captures.

## Pitfalls (each cost real time; do not relearn)

- **Thermal σ_nγ/σ_np = 1.0×10⁻⁸ is sub-keV ONLY** — it rises to ~10⁻⁴ at
  MeV. Using it everywhere underreports MeV rates ×300 (this bug shipped in
  a first pass of the MeV analysis and was caught on re-audit). Use direct
  (n,γ) counts.
- **n/pulse anchors:** sub-keV runs 7.31×10⁶; fullrange runs 2.2628×10⁷
  (integral of `flux_n_pulse_NOisolet_100bpd` in
  `data/fluxEAR2-Ph3_in_different_units.root`). Alberto's table assumes
  3.29×10⁷ (31% hotter) — compare per incident neutron, never per pulse.
- **³He(n,p)t = `neutronInelastic`**, not `nCapture`; `nCapture` in He3Gas
  is the X17 channel. Wall volume = `He3Cap_Al`, CFRP = `He3Cap_CFRP`.
- **X17 branching 2.5% is per IPC pair** (× α_IPC = 2.1e-3 per capture),
  not per capture.
- Numerical accident: α_IPC × BR_X17 × pulses/day = 1.01 → capt/pulse ≈
  X17/day. Handy for sanity checks.
- uproot string branches: `np.asarray(..., dtype=object)` immediately.
- lxplus: `ssh lxplus` lands in `$HOME`; repo at
  `/afs/cern.ch/work/d/dneff/git/MX17_Full_Geant`; env via
  `source scripts/setup_lxplus.sh`. Local python:
  `source ../../PycharmProjects/nTof_x17/.venv/bin/activate`.

## Caveats that bound the final numbers

- **At-rest kinematics:** pairs were generated with E* = 20.58 MeV, no CM
  boost. MeV captures have E* ≈ 20.58 + ¾E_n (+1.5 MeV at E_n = 2 MeV) —
  opening angle, energy sharing, TOF all shift. The pools are the right
  machinery baseline; the **generator E_n-kinematics extension + new run A**
  are needed before final acceptance/sensitivity numbers.
- α_IPC = 2.1×10⁻³ and X17/IPC = 2.5% are unconfirmed assumptions
  (Alberto follow-up open). Everything scales linearly in both.
- Above 20 MeV G4NDL ends — that decade is model-dependent (1% of budget).

## Open items beyond Stage 4

- Generator extension for E_n-dependent kinematics → re-plan run A.
- Run C (biased wall-γ source) for calorimeter spectra of wall backgrounds;
  run D (single-particle trigger curves).
- Alberto: α_IPC provenance, X17/IPC basis, flux-normalisation
  reconciliation, capsule-geometry update of the table.
- γ-flash recovery measurement at EAR2 (experimental, blocks the
  flash-parameter choice in the sampler).
