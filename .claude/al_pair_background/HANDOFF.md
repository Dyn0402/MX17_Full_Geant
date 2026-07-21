# Handoff: Aluminium pair-production background to the X17 / IPC measurement

**Created:** 2026-07-21 · **Priority: HIGH — potential show-stopper for the thermal X17/IPC measurement**
**Companion docs:** [../thermal/thermal_campaign_handoff.md], [../../CAMPAIGN_STATUS.md],
[../e0_branch/ipc_roadmap.md] (IPC multipolarity), [../../analysis/thermal_2cm/] (this session's figures)

---

## TL;DR

The **²⁷Al(n,γ)²⁸Al capture γ (7.72 MeV) pair-produces** in the Al capsule, making a
**real e⁺e⁻ pair — the same final state as our IPC and X17 signal.** In the thermal
gate this is **~31 pair-legs/pulse ≈ 6×10⁵/day** (trigger level) against an IPC signal
of **~1.4/day** and an X17 signal of **~0.035/day** — a raw ratio of order **10⁵–10⁶ : 1**.

The two rejection handles we originally assumed are **both unavailable**:
- **Vertex** — resolution is insufficient to reject a pair originating anywhere near the
  capsule (the capsule sits at r ≈ 1–2.6 cm, right around the r < 1 cm gas signal region).
- **Energy** — the LS calorimetry is not working, so we cannot cleanly separate the
  7.72 MeV Al pair from the 16.8–20.58 MeV signal by total energy.

**The only surviving discriminant is the e⁺e⁻ opening angle from the MM tracking**
(plus, possibly, track-level momentum handles — see Task 4). This handoff specifies a
thorough analysis to (a) characterize the Al pair-production background fully, (b)
compute its opening-angle distribution vs IPC and X17 **assuming infinite statistics**,
and (c) determine whether the thermal X17/IPC measurement survives it.

---

## Why this is serious (what changed)

Earlier this session we reassured ourselves that the Al pair background is discriminable
by **energy** (7.72 ≪ 16.8–20.58 MeV) and **vertex** (capsule vs gas). Both of those
assumptions are now off the table given the real detector (no working calorimetry, no
vertex resolution near the capsule). That removes the clean handles and promotes this
from "a trigger nuisance we reject offline" to **"a physics background in the same
observable space as the signal, at ~10⁵× the rate."**

The pair is kinematically distinct from the signal only in the **opening angle** and in
the **individual lepton energies** (Al leptons are few-MeV; signal leptons up to ~19 MeV).
Whether those survive our resolution — and whether the enormous rate leaks into the signal
opening-angle window through tails / multiple scattering / mismeasurement — is the whole
question.

---

## What we established this session (anchor numbers)

Geometry: current `_2cm` (backscint 2.0 cm). In-gate = thermal gate, t > 1 ms ⇔ E_n < 2 eV.
Trigger "leg" = SiPM-wall AND plastic ≥ 0.5 MIP in one arm. MIP: SiPM 458 keV, plastic
≈ 3.47 MeV (2.0 cm, scaled from the 2.5 cm muon MPV 4.334 MeV — **needs a proper 2.0 cm
muon calibration**).

- **Al legs in-gate: 197 / pulse** (of 205 total legs; Al = 96% of all legs).
- **Pair-production fraction of Al legs: 15.6%** → **~31 pair-legs / pulse ≈ 6×10⁵ / day.**
  (84.4% of Al legs are double-Compton; classifier = positron energy in the deposits,
  `scripts/leg_mechanism_split.py`. Median e⁺ energy fraction in pair legs = 59%.)
- Background energy deposition is **~99.5% e⁻/e⁺** (photon-induced), <0.3% proton recoil
  (`scripts/count_particle_composition.py`).
- The pair-production **vertex is mostly in the Al capsule** (`conv` in `He3Cap_Al`);
  the 7.72 MeV γ converts, then the e⁺e⁻ shower deposits in the detectors
  (real event: `analysis/thermal_2cm/event_display_al_leg.png`, eventID 126734).
- **Signal for comparison:** IPC in-gate **6.58×10⁻⁵ /pulse = 1.39 /day** (2cm reweight,
  `analysis/reweight_2cm/`); X17 = 2.5% of IPC ≈ 0.035/day.
- **Raw ratio Al-pair : IPC ≈ 4×10⁵ : 1** (trigger level; the MM-acceptance ratio is a
  key deliverable, see Task 3 — it will differ but stay very large).

> **Caveat on "31/pulse":** this is the SiPM∧plastic *leg trigger* rate. The signal is
> defined by **MM pair reconstruction** (two tracks from a vertex). The physics-relevant
> number is the rate of Al pairs **reconstructed as a pair in the MM acceptance**, which
> the analysis must compute directly (Task 3), not the leg rate.

---

## The physics of the three opening-angle distributions (what to expect)

All three are e⁺e⁻ final states; they differ in origin and kinematics:

| source | origin | total energy | opening-angle shape (expected) |
|---|---|---|---|
| **Al pair** | external γ→e⁺e⁻ conversion, 7.72 MeV, in the Al capsule | ~7.72 MeV | **forward-peaked, small angle** (θ ~ m_e c²/E ~ few°), **but** with a large-angle tail from (a) asymmetric-energy pairs and (b) heavy multiple scattering of the soft leptons |
| **IPC** | internal pair conversion of the 20.58 MeV ³He(n,γ)⁴He transition | ~20.58 MeV | **monotonically falling** from small angle; shape depends on **multipolarity (E0/M1/E2)** — coordinate with [../e0_branch/] |
| **X17** | ⁴He\* → ⁴He + X17(16.8 MeV), X17 → e⁺e⁻ | ~20.58 MeV | **peak at a large, characteristic angle** (ATOMKI-like); position set by the X17 mass and the 20.58 MeV kinematics |

**The central tension:** X17 lives at *large* opening angle, Al at *small* — so they
separate *in principle*. But Al outnumbers the signal by ~10⁵, so **even a 10⁻⁴–10⁻⁵
large-angle tail of the Al distribution equals the entire X17 signal.** The two things
that fill that tail and MUST be quantified with detector realism:
1. **Multiple Coulomb scattering** of the few-MeV Al leptons in the capsule + MM +
   walls — softer leptons scatter far more than the ~19 MeV signal leptons, smearing
   the reconstructed opening angle to larger values.
2. **Mismeasurement / wrong-track pairing** — combinatorial pairing of tracks from
   different particles (huge singles rate: ~2000 SiPM / 900 plastic singles/pulse).

---

## Constraints (given by the experiment — assume these are firm)

- **No vertex-based rejection** for anything near the capsule (r ≲ 3 cm). Do NOT rely on
  "pair originates in the gas" as a cut. (You may still study the *radial* distribution
  as a weak statistical handle — capsule is a shell at r≈1–2.6 cm, signal a volume at
  r<1 cm — but assume it cannot be used as a hard cut.)
- **No calorimetry** (LS essentially absent). Do NOT assume a total-energy measurement
  from scintillator light. Any energy information must come from the **MM tracking
  itself** (multiple-scattering-inferred momentum, range) — see Task 4.
- Thermal gate only (t > 1 ms). Epithermal/prompt is a separate (larger) regime.

---

## The analysis — tasks

**Task 0 — reproduce/verify the anchor numbers.** Rerun `leg_mechanism_split.py` and the
capture/particle tallies on `neutrons_thermal_trig_2cm` to confirm the 197 Al legs/pulse,
16% pair fraction, etc. Get a proper **2.0 cm plastic MIP calibration** from a muon run
(current value is scaled, ±20% on the plastic threshold → affects rates).

**Task 1 — extract the Al pair-production truth kinematics.** The neutron sim does NOT
save the e⁺e⁻ truth for conversion pairs. Two options:
  (a) **From the trajdump** — at each `conv` vertex in `He3Cap_Al` (and other volumes),
      read the birth momenta/energies of the daughter e⁺ and e⁻ → opening angle, individual
      KE, total energy, vertex radius. Tools staged: `scripts/find_al_leg_event.py` shows
      the parsing; a 150k-event trajdump sits at `lxplus:/tmp/dneff/trajbig_traj_t0.csv`
      (296 MB, **ephemeral — regenerate with `--trajdump` if gone**). This is the fastest
      path to a first distribution.
  (b) **Add truth branches to the sim** — cleaner and higher-stats: in `SteppingAction`/
      `EventAction`, when a `conv` produces an e⁺e⁻ in a capture context, record the pair
      (like `PairKinematics` for the signal). Then run a dedicated (optionally
      nCapture-biased, see `--bias-ncapture`) campaign for high statistics. Recommended
      for the final result.
  Extend beyond Al to **all pair-producing capture γs**: N(n,γ) **10.8 MeV** (closest to
  signal!), Si 8.47, Cu 7.9, C 4.95, ³He(n,γ) 20.58 (this is the signal). Rank by rate ×
  proximity to the signal opening-angle window.

**Task 2 — get the IPC and X17 truth kinematics.** Directly from `pairs_thermal_trig_2cm`
(EventTree branches `openingAngle`, `em_ke`, `ep_ke`, `inv_mass`, `vtx_*`, `event_type`
0=X17 1=IPC). No new sim needed. **Confirm the IPC multipolarity assumption** in
`X17PrimaryGenerator` against [../e0_branch/] — the IPC opening-angle shape depends on it.

**Task 3 — rate normalization in the MM acceptance.** Compute, per pulse and per day:
Al pairs (and other-material pairs) **reconstructed as a 2-track pair in the MM
acceptance**, vs IPC and X17. This replaces the trigger-level 31/pulse with the
analysis-level number and gives the true S/B before any opening-angle cut.

**Task 4 — the opening-angle distributions (the headline deliverable).**
  (a) **Infinite-stats truth shapes** (as the user asked): overlay normalized
      opening-angle distributions for Al, IPC, X17. Shows intrinsic separability.
  (b) **Rate-weighted**: multiply each by its Task-3 rate → the distribution we'd actually
      measure (Al will dominate; see where the X17 peak sits relative to the Al tail).
  (c) **With detector realism**: apply MM acceptance + angular resolution + **multiple
      scattering** (this is the crux — propagate the leptons, or fold a realistic
      per-track scattering model that depends on lepton energy). The Al leptons (few MeV)
      scatter far more than signal leptons (~19 MeV); quantify how far their tail reaches
      into the X17 window.
  Also explore whether the **individual lepton energies / track quality** (inferred from
  MM multiple scattering ∝ 1/pβ, or range) give an additional handle to suppress the soft
  Al leptons without a calorimeter.

**Task 5 — verdict.** Given no vertex, no calorimetry: what opening-angle (± track-energy)
selection gives what residual Al background under the X17 peak and under the IPC spectrum?
Is the thermal X17/IPC measurement feasible, and at what significance? Feed the residual
into the significance projection (see `docs/slides/`, [../thermal/]).

---

## Data & code pointers

- **Neutron sim (Al captures + pair production):**
  `/eos/experiment/ntof/data/x17/full_sim/neutrons_thermal_trig_2cm` (10⁹, validated
  100/100). Also `neutrons_epi_trig_2cm` (5×10⁸) for the t<1 ms regime.
- **Signal pairs (IPC + X17, truth kinematics):**
  `/eos/experiment/ntof/data/x17/full_sim/pairs_thermal_trig_2cm` (10⁷, 50/50 X17/IPC,
  vertices from the measured thermal self-shielding profile).
- **Build:** `source scripts/setup_lxplus.sh && bash scripts/build.sh` (Geant4 11.2 /
  ROOT 6.32 from CVMFS). Executable supports `--trajdump N` and `--bias-ncapture FACTOR`.
- **Session scripts (all committed):** `scripts/leg_mechanism_split.py` (pair vs Compton
  classifier), `count_particle_composition.py`, `find_al_leg_event.py`,
  `plot_event_display.py`, `count_timedist_*thermal.py`, `plot_reaction_table.py`.
- **Session figures:** `analysis/thermal_2cm/` — `event_display_al_leg`,
  `leg_mechanism_schematic` (84%/16%), `leg_mechanism_split`, `reaction_table`,
  `particle_composition`, `timedist_thermal`, `timedist_bysource`.
- **Trigger analysis / MIP calibration:** `scripts/analyze_trigger_thermal.py`,
  `analysis/trigger_thermal/trigger_scan.json` (2.5 cm; re-run on `_2cm` for consistency).

---

## Pitfalls & subtleties

- **Leg rate ≠ pair rate.** 31 pair-legs/pulse is the SiPM∧plastic trigger rate. Compute
  the MM-reconstructed pair rate for the physics comparison (Task 3).
- **Multiple scattering is the whole ballgame** for the tail. The Al leptons are soft
  (few MeV); a truth-only opening angle will look cleanly separated and give a false sense
  of security. The realistic (scattered) distribution is what determines feasibility.
- **Plastic MIP for 2.0 cm is scaled, not measured** (±20%). Affects all rates.
- **IPC multipolarity** sets the IPC opening-angle shape — confirm the generator
  assumption (E0/M1/E2) with the E0 branch work before quoting IPC/Al separation.
- **Other pair backgrounds:** N(n,γ) at 10.8 MeV gives higher-energy pairs (closer to the
  signal window) — do not ignore it just because Al dominates by rate.
- **Combinatorial pairing:** with ~2000 SiPM + ~900 plastic singles/pulse, wrong-track
  pairing can fake large opening angles independent of the true pair kinematics.
- **Don't reuse the June (pre-overhaul) geometry** for anything — use `_2cm`.
- **Neutron HP transport:** ³He(n,p)t is `neutronInelastic`; the radiative/γ-emitting
  captures are `nCapture` (Al, H, etc.). The signal ³He(n,γ) is also `nCapture` in
  `He3Gas`.

---

## Deliverables

1. Al pair-production **rate in the MM acceptance** (per pulse, per day), and the same for
   IPC and X17 → the true S/B before cuts.
2. **Opening-angle distributions**: (a) infinite-stats truth shapes (Al/IPC/X17 overlaid,
   normalized), (b) rate-weighted, (c) with MM resolution + multiple scattering.
3. **Residual Al background** under the X17 peak and IPC spectrum after the best
   opening-angle (± track-energy) selection.
4. **Feasibility verdict** for the thermal X17/IPC measurement given no vertex and no
   calorimetry, folded into the significance projection.
5. Extension to the **other pair-producing capture γs** (esp. N 10.8 MeV).
