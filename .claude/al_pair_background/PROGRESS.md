# Al pair-production background — analysis progress

**Started:** 2026-07-21 (this session, continuing the HANDOFF.md task). Author: Claude (Opus 4.8) + Dylan.

## What's done

### Sim: γ→e⁺e⁻ conversion truth branch (Task 1b — the recommended high-stats path)
Added a compact `ConvPairTree` to the simulation that records the birth truth of
**every** photon `conv` (pair production) during neutron transport:
- `include/EventData.hh` — new `ConvPair` struct + `convPairs` vector on `EventData`.
- `src/SteppingAction.cc` — `RecordConvPair()`: on a gamma `conv` step, reads the
  e⁺/e⁻ secondaries (`GetSecondaryInCurrentStep`) and stores birth KE, momentum
  unit vectors, opening angle, conversion vertex + volume, and γ energy.
- `src/RunAction.cc` — `ConvPairTree` (flat, one entry per conversion): `eventID,
  weight, neutron_E_eV, capture_vol, gamma_E, conv_vol, vx/vy/vz, em_ke+dir,
  ep_ke+dir, openingAngle`.
- `include/SteppingAction.hh` — declaration.

Built on lxplus (`/afs/cern.ch/work/d/dneff/git/MX17_Full_Geant`), validated on a
20k test: `gamma_E = 7.725 MeV` exactly = the ²⁷Al(n,γ) 7.72 line, `conv_vol =
He3Cap_Al`. **Local repo source files are edited but NOT yet committed.**
(NB: `-t 4` MT aborts on the lxplus *login* node — a cling+MT quirk; batch jobs
run single-thread and are fine.)

Also added `gamma_trackID` to `ConvPairTree` (the converting γ's trackID): the
daughter leptons' MM hits are the `HitTree` DriftGas hits with
`parentID == gamma_trackID`, split e-/e+ by the `particle` field — this is the
robust join for the reconstructed-opening-angle analysis (daughter trackIDs
aren't yet assigned inside the conv step, so we anchor on the γ). Validated.

### Campaign (running — resubmitted with gamma_trackID)
`neutrons_convpair_2cm` — 100 jobs × 10⁷ thermal neutrons (E=[0.001,2] eV), the
**exact config of the validated `neutrons_thermal_trig_2cm`**, writing to
`/eos/experiment/ntof/data/x17/full_sim/neutrons_convpair_2cm/`.
Jobdir `/afs/cern.ch/user/d/dneff/condor/mx17_convpair_2cm`, seed 4217, cluster
11817797. Analog Al conv rate ≈ 7×10⁻⁵/event (thermal window) → expect ~7×10⁴ Al
conv pairs from 10⁹; ~half reach the MM, ~20% form a 2-track MM pair.
**Status: running (first submission was killed+resubmitted to add gamma_trackID;
no data lost — 0 jobs had finished).** Completion check: `grep -l "Job done"
<jobdir>/logs/*.out` (NOT file size — ROOT autoflushes HitTree mid-run, so
>0.5 GB does NOT mean complete; ConvPairTree only lands at end-of-job).

### Signal truth (Task 2 — done, no new sim)
From `pairs_thermal_trig_2cm` EventTree (`scripts/extract_signal_openingangle.py`):
- **X17 opening angle: median 117°, IQR [111,130]** — sharp kinematic threshold
  at ~108° (m=16.8 boson from 20.58 MeV), falling to 180°.
- **IPC opening angle: median 30°, IQR [11,71]** — falls from small angle, but a
  **14.3% tail sits above 100°** (under the X17 peak). This is the IPC-vs-X17
  separation problem, independent of Al.
- Figure: `analysis/al_pair/signal_openingangle.png` (+ npz).

## Key findings so far

1. **Soft-lepton energy is a powerful, nearly-orthogonal discriminant.** The Al
   pair's **softer lepton is always ≲3 MeV**; the X17 signal's **softer lepton is
   always ≳4 MeV** (it is massive → constrained energy sharing). In truth these
   **do not overlap.** If the MM can infer even a coarse lepton momentum (from
   multiple-scattering ∝ 1/pβ or range), Al rejection vs X17 could be far cleaner
   than opening angle alone. IPC's softer lepton spans 0–10 MeV, so this handle
   helps X17 more than IPC. (Panel B of the figure.)

2. **The generator's IPC has NO multipolarity.** `X17PrimaryGenerator::GenerateIPC`
   samples `Mee ∝ 1/Mee` on [2mₑ, 20.58] then uses the **same isotropic two-body
   decay + boost** as X17 — no E0/M1/E2 angular correlation. So the IPC
   opening-angle curve (and its 14% >100° tail) is a **modelling baseline, not a
   first-principles IPC distribution**. This is the "IPC shape systematic" that
   dominates the real CL — must be reconciled with the E0 branch before quoting
   any Al/IPC separation. See `.claude/e0_branch/`.

3. Trajdump fast-path gave only ~22 Al pairs (median θ 27°, inflated by first-step
   MSC) — too few and biased; confirms the ConvPairTree campaign is the right call.

4. **MSC smears the Al opening angle up by tens of degrees (preview, n=3).** The
   reconstruction pipeline (`scripts/analyze_convpair_reco.py`, MM-entry direction
   = momentum dir of earliest DriftGas hit) is validated on a 200k test: Al truth
   θ median ~15° but **reconstructed θ median ~60°** (median |Δ|~52°). Direction
   and magnitude match expectation (few-MeV e± through the Al capsule). Also ~⅔ of
   MM pairs enter the *same arm* → 2-track separation is a real challenge. Needs
   full-campaign stats to quantify the tail past 108° (X17 threshold).

## REAL results (20 complete files, 200M events, 30,165 conv pairs) — 2026-07-21

Rate weight = 4 292 400.81 / 2e8 = 0.02146 /pulse; ×1.929e4 pulses/day.

**Al capsule conv (He3Cap_Al): 14,371 pairs = 308/pulse = 5.95×10⁶/day produced.**
- Truth opening angle: **median 22°**, **2.37% above 108° (X17 threshold)** →
  1.4×10⁵/day in the X17 angular region *at truth, produced level* (vs X17
  0.035/day). γ line 7.73 MeV confirmed (+ cascade lines 1.8/3.0/4.1/4.3/4.7).
- **Softer-lepton KE: median 1.31 MeV, only 0.007% above 4 MeV** — X17's softer
  lepton is a 4–10 MeV box (massive boson), so a ~4 MeV soft-lepton cut removes
  ~99.99% of Al with ~zero X17 loss, IF the MM can measure lepton momentum.
- Vertex radius median 2.7 mm [1.1, 9.4] — capsule shell, overlaps r<10mm signal.

**The Al 7.72 MeV γ converts everywhere, not just the capsule:** counting by the
7.7 MeV line (all conv volumes) → 15,342 pairs = 6.35×10⁶/day (He3Cap_Al 8087 +
LiqScint 2033 + walls…). Conv location matters for MSC (capsule = longest
scattering path before the MM).

**Per-material ranking (produced /day):** Al capsule 5.95e6, LiqScint_1 1.97e6,
LS_VesselCFRP 7.7e5, BackScint L+R 1.2e6, He3Cap_CFRP 6.0e5, PlasticScint 3.9e5,
PCBs/Micromesh/cathode ~1e5 each.

**N(n,γ) 10.8 MeV (handoff flag):** 58 pairs in [9.5,11.5] MeV = 2.4×10⁴/day,
mostly LiqScint/CFRP nitrogen. 250× below Al and MORE forward (median 9°, harder
leptons weaken the soft-lepton handle) — subdominant but real. Si ~8.5 line
8.3e3/day.

**RECO (MSC) preview from 200k test (n=3):** truth θ 15° → reco θ 60°. Full-stats
reco analysis running (`analyze_convpair_reco.py`) — the decisive tail number.

Figures: `analysis/al_pair/convpair_truth.png` (3-panel: opening angle, soft-lepton
KE, rate-weighted), `signal_openingangle.png`.

## Reconstructed (MSC) opening angle — Task 4c/3, ONE MM definition for all species

MM 2-track "pair" = both leptons make ≥1 DriftGas hit; entry direction =
momentum dir of earliest-time DriftGas hit (captures scattering through the Al
capsule + air before the tracker). `analyze_convpair_reco.py` (Al, parentID==
gamma_trackID), `analyze_signal_reco.py` (X17/IPC primaries, parentID==0).
NB run reco on **condor** (batch) — the login node kills the ~22-min job; it
writes npz to EOS `convpair_out/`. Signal reco (10 pairs files, 300k events)
ran fine as a login nohup.

**Signal reco (300k events):**
- X17: MM acc **33.0%**; truth θ med 117° → reco 120° (hard leptons, little MSC);
  reco **>108° = 67%** (33% scattered below threshold); same-arm only 4%.
- IPC: MM acc **32.0%**; truth θ med 24° → reco 42°; reco **>108° = 14%**;
  same-arm 55%.

**Al reco (single-file preview n=715, full condor run pending):** MM acc ~9.8%;
truth θ med 17° → reco 39°; reco **>108° ≈ 10%**; same-arm 63%.

**Preliminary rate picture in the X17 window (reco θ>108°), per day:**
- Al ≈ 5.95e6 × 0.098 × 0.10 ≈ **6×10⁴/day**
- IPC ≈ 1.39 × 0.32 × 0.14 ≈ 0.06/day
- X17 ≈ 0.035 × 0.33 × 0.67 ≈ 0.008/day
→ **Al : X17 ≈ 10⁷ : 1 from opening angle ALONE — hopeless.** The soft-lepton
handle (both leptons >4 MeV; Al 0.007% survive) is essential: it brings Al to
~6/day, but that is still ~800× the X17 rate (and needs a real MM momentum
resolution, treated as ideal here). **Leaning: the thermal X17 measurement is not
feasible against the Al pair background with opening angle + ideal soft-lepton
alone; IPC (0.06/day in window) is far more favourable.** Confirm with full Al
reco stats + momentum-resolution realism + the IPC-multipolarity caveat before
finalizing (Task 5). Verdict script: `scripts/plot_convpair_verdict.py`.

## Next steps (when campaign data lands)
- Run `scripts/analyze_convpair_truth.py` → Al truth opening angle / softer-lepton
  KE / vertex-r, per-material table, **rate/pulse & /day** (weight/pulse =
  4 292 400.81 / N_sim; ×1.929e4 pulses/day). Rank materials incl. **N 10.8 MeV**.
- **Task 3 / 4c (needs HitTree):** reconstruct the e⁺e⁻ opening angle from the MM
  `DriftGas` hits in conv-pair events → the REAL detector-level (scattered)
  distribution — no analytic MSC model needed, the sim already transported the
  leptons through the true geometry. Compare truth vs reconstructed; measure the
  Al large-angle tail into the X17 window; MM-acceptance rate normalization.
- **Task 5 verdict:** best opening-angle (± soft-lepton/track-momentum) selection,
  residual Al under the X17 peak and IPC spectrum, feasibility + significance.

## Data & code
- Sim source (edited, uncommitted locally; scp'd + built on lxplus).
- `scripts/extract_signal_openingangle.py` (Task 2 + trajdump cross-check)
- `scripts/analyze_convpair_truth.py` (Task 1b/3/4a truth — staged, run when data lands)
- `scripts/plot_signal_openingangle.py` (local figure)
- `analysis/al_pair/` — figures + npz (local)
