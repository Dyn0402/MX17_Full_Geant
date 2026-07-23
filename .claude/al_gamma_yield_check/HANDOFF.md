# Handoff — Al(n,γ) 7.72 MeV gamma yield: analytic cross-check vs Geant4

**Written 2026-07-23 by Claude (Opus 4.8) for Dylan. This is a task specification, not a result — no calculation has been done here.**

---

## 1. The problem

We have a Geant4 number for the ²⁷Al(n,γ) 7.72 MeV capture-gamma rate from the
He-3 capsule's aluminium pressure vessel. Separately, a back-of-envelope
calculation was done using **13 g of aluminium and the same EAR2 beam from 1 eV
down**, and **the gamma yield it gave looked low**.

That hand calculation used (or is suspected to have used) a **thin-disk /
single-pass approximation**: treat the Al as a slab of some areal density
`N·t [atoms/cm²]`, take capture probability `P = N·t·σ_nγ`, multiply by the
in-window neutron flux. The worry is that this approximation is not appropriate
for our actual geometry.

**Your job:**

1. Redo the analytic estimate properly, and **compare thin-disk against a less
   approximate treatment** — quantify how much the thin-disk assumption costs
   and in which direction.
2. Produce a **simplified but defensible calculation** of the expected Al(n,γ)
   rate that we can put in a note.
3. **Compare it to what Geant4 actually gives** and explain any discrepancy.

Do not just reproduce the low number — the point is to find out whether the low
number is right, and if not, why.

---

## 2. Why thin-disk is suspect here (the things to check)

These are hypotheses, not conclusions. Test the ones that matter.

- **The vessel is not a slab of uniform thickness.** Al thickness along the beam
  varies by an order of magnitude with radius: **5.5 mm on-axis at the nose**
  vs **0.6 mm through the side wall**. A single effective thickness is a fiction;
  the right thing is a thickness-vs-radius profile folded against the beam's
  radial profile. Both are given below.
- **Which 13 g the beam actually sees.** The full 13 g includes the neck/valve
  stem at z = +40…+51 mm, which sits *behind* 17 cm³ of 500 bar He-3. He-3 is
  enormously absorbing at thermal energies (σ_np ≈ 5333 b at 0.0253 eV, 1/v) —
  the gas is optically thick below ~1 keV (see `docs/he3_self_shielding_note.md`,
  optical depth ~150). So a large fraction of the mass may be **shadowed** and
  contribute far less than its mass share. This pushes the estimate *down*, not
  up — but it must be accounted for before any comparison is meaningful.
- **1/v and the sub-thermal flux.** The beam runs down to 1 meV. σ_nγ for ²⁷Al is
  1/v below the resonances, so at 1 meV it is ~5× the 0.0253 eV value. **If the
  hand calc used σ(0.0253 eV) with the whole sub-eV flux, it underestimates.**
  This is a prime suspect for "looked low." Do the integral over the flux
  spectrum, not at a single energy.
- **Elastic scattering lengthens the path.** ²⁷Al elastic σ ≈ 1.4 b vs capture
  ≈ 0.23 b at thermal — scattering dominates by ~6:1. Neutrons random-walk inside
  the Al and re-enter it after scattering, so the effective path exceeds the
  geometric one. Single-pass ignores this entirely.
- **Albedo from the surrounding detector.** The capsule sits inside plastics, a
  liquid-scintillator layer, and MM structure — all hydrogenous moderators that
  scatter thermal neutrons back onto the capsule. Geant4 sees this; a single-pass
  slab calc does not. Also a candidate for the shortfall.
- **~12 % of the thermal beam misses the capsule entirely** (see §4) and can only
  reach the Al via albedo.

The honest expectation is that some of these push up and some push down. Rank
them by size rather than listing them.

---

## 3. Geometry — the Al vessel as modeled

All from the `He3Cap_Al` and `He3Gas` polycones in
`src/DetectorConstruction.cc:320-421`. STEP source: "MASTINU X17 HPRV 00 01".
Polycone axis is local z, rotated −90° about x, so **the capsule axis is the beam
axis (+Y in world)**. z below = distance along the beam, z<0 = upstream.

| Quantity | Value |
|---|---|
| **Al mass** | **13.24 g** (verified: 4.9068 cm³ × 2.699 g/cm³, G4_Al) |
| Al volume | 4.9068 cm³ (= 21.9095 cm³ outer polycone − 17.0027 cm³ gas cavity) |
| Total length | 86 mm (z = −35 mm tip → +51 mm valve top) |
| Barrel OD | 21.2 mm (r = 10.6 mm), flat over z = −21…+21 mm |
| Neck / valve OD | 7.0 mm (r = 3.5 mm), z = +40…+51 mm |
| **Gas bore** | **r = 10 mm, L = 40 mm** (z = −20…+20 mm) |
| Gas cavity full extent | z = −29.5 → +50.7 mm; hemispherical lower dome; necks to r = 0.75 mm fill channel |
| **Barrel wall thickness** | **0.6 mm Al** (r = 10.0 → 10.6 mm) |
| **Nose thickness on-axis** | **5.5 mm Al** (z = −35 → −29.5 mm) |
| CFRP wrap | 0.9 mm over the Al outer surface (`He3Cap_CFRP`, separate volume — not Al, but it moderates and is itself a capture site) |
| He-3 fill | 17.0 cm³ @ 500 bar, ρ = 62.7 mg/cm³ → ≈ 1.07 g |

Full polycone vertex arrays (z, r_outer) for Al, gas, and CFRP are in
`src/DetectorConstruction.cc:341-412` — **copy them from there rather than
re-typing**, and use them to build the thickness-vs-radius profile. A mirror of
the same profile lives in `scripts/plot_geometry.py`.

²⁷Al is 100 % naturally abundant, A = 26.9815 g/mol → the atom count follows
directly from 13.24 g.

---

## 4. Beam — normalization and radial profile

The gun samples energy from `flux_n_pulse_NOisolet_100bpd` in
`data/fluxEAR2-Ph3_in_different_units.root` and transverse position from
`Lambda2D` in `data/lamda2DvsEn_EAR2.root` (240 log-E rows × 3000 radial bins to
r = 3 cm). Implementation: `src/X17PrimaryGenerator.cc:355-383`. Beam is along
+Y; row weights already include the 2πr area factor, so the radial CDF is
sampled directly in r.

**Flux normalization (integrals of the flux histogram, verified):**

| Window | n/pulse |
|---|---|
| Full range | 2.263×10⁷ |
| < 1 keV | 7.312×10⁶ |
| **< 2 eV** (the >1 ms thermal gate) | **4.283×10⁶** |
| **< 1 eV** (the window the suspect calc used) | **3.996×10⁶** |
| < 0.5 eV | 3.698×10⁶ |

**Pulses/day = 1.929×10⁴** (7×10¹² ppp). These anchors are consistent with
`CAMPAIGN_STATUS.md:17` and `.claude/mev/HANDOFF_MEV_ANALYSIS.md:75-80`.

**Radial profile, flux-weighted (verified from `Lambda2D`):**

| Window | r₅₀ | r₆₈ | r₉₀ | r₉₅ |
|---|---|---|---|---|
| Thermal (<0.5 eV) | 0.26 cm | 0.47 cm | 1.33 cm | 1.86 cm |
| < 1 keV | 0.21 cm | 0.36 cm | 1.10 cm | 1.63 cm |

Fractions by radial zone (thermal <0.5 eV / <1 keV):

| Zone | thermal | <1 keV |
|---|---|---|
| Through the gas bore (r < 1.00 cm) | 85.0 % | 88.7 % |
| Al barrel-wall annulus (1.00–1.06 cm) | 1.1 % | 0.9 % |
| CFRP annulus (1.06–1.15 cm) | 1.5 % | 1.1 % |
| Misses the capsule (r > 1.15 cm) | 12.4 % | 9.3 % |

**Two consequences for the calculation.** First, only ~1 % of the beam hits the
barrel side wall — the Al the beam actually traverses is overwhelmingly the
**5.5 mm nose**, not the 0.6 mm wall. Any thin-disk calc using the wall thickness
is wrong by ~10×. Second, `Lambda2D` is **truncated at r = 3 cm**, so the "misses"
fraction is relative to what is in the histogram; halo beyond 3 cm is not modeled
in the sim at all, and is therefore also absent from the Geant4 number you are
comparing against. Do not add it to the analytic side only.

---

## 5. Nuclear data you need

**Verify all of these against ENDF/B before using — they are quoted here from
memory as starting points, not as vetted inputs.**

- ²⁷Al(n,γ) thermal capture cross section ≈ **0.231 b at 0.0253 eV**, 1/v below
  the first resonance (5.9 keV). Extend as `σ(E) = σ_th·√(0.0253/E)` over the
  sub-eV window.
- ²⁷Al elastic ≈ 1.4 b thermal (needed for the path-lengthening estimate).
- The **7.724 MeV line is the ground-state transition of ²⁸Al** and is only a
  *fraction* of captures — the cascade also produces 1.8/3.0/4.1/4.3/4.7 MeV lines
  (all five observed in our sim, `.claude/al_pair_background/PROGRESS.md:88`).
  **Get the per-capture intensity of the 7.724 line from a capture-gamma library
  (IAEA/Lone) — do not assume 1 γ per capture.** Getting this wrong is another
  strong candidate for a factor-level discrepancy.
- ³He(n,p)t ≈ **5333 b at 0.0253 eV**, 1/v — needed for the gas-shadowing term.
  Tabulated in `data/He3.h5` (OpenMC, ENDF); `scripts/reweight_ipc_vs_time.py:49`
  shows how to read a reaction MT out of it.

---

## 6. What Geant4 says (the comparison target)

From the full-statistics conv-pair campaign, documented in
`.claude/al_pair_background/PROGRESS.md:83-99` and `VERDICT.md`:

- Rate weight: **4 292 400.81 / 2×10⁸ = 0.02146 per pulse**, × 1.929×10⁴ pulses/day.
- **Al capsule γ→e⁺e⁻ conversions (`He3Cap_Al`): 14,371 pairs = 308/pulse = 5.95×10⁶/day produced.**
- Counting by the 7.7 MeV line across *all* conversion volumes: 15,342 pairs =
  6.35×10⁶/day (He3Cap_Al 8087 + elsewhere).
- Confirmed line energy: `gamma_E = 7.725 MeV`.

**Careful — that is a *pair-conversion* rate, not a γ-production rate.** It is
already folded with the probability that the γ pair-produces in a given volume.
To compare against an analytic γ yield you need one of:

- **(preferred) Count captures directly.** `EventTree` carries `capture_vol` and
  `capture_proc` per event (`src/RunAction.cc:48-49,138-139`). Selecting
  `capture_vol == "He3Cap_Al" && capture_proc == "nCapture"` gives the **number of
  Al radiative captures** with no pair-production factor in the way. This is the
  clean apples-to-apples quantity — **use this.**
- Or unfold the conversion probability from the pair number, which is messier and
  geometry-dependent.

**Data locations (lxplus/EOS):**
- Thermal campaign: `/eos/experiment/ntof/data/x17/full_sim/neutrons_thermal_trig_2cm` (10⁹ neutrons, 1 meV–2 eV)
- Conv-pair campaign: `/eos/experiment/ntof/data/x17/full_sim/neutrons_convpair_2cm`

Example readers to copy: `scripts/count_timedist_thermal.py`,
`scripts/analyze_convpair_truth.py`. Environment: `source scripts/setup_lxplus.sh`
on lxplus, or `source ../../PycharmProjects/nTof_x17/.venv/bin/activate` locally.

One useful lever already in the sim: `SimConfig::disableAlCapsule` swaps the Al
vessel for vacuum (`include/SimConfig.hh`, `src/DetectorConstruction.cc:393-395`),
and a no-Al run already exists at `analysis/al_pair_crosscheck/timedist_noAl.npz`.
Differencing Al vs no-Al isolates the Al contribution if you need it.

---

## 7. Suggested deliverables

Write results into this directory (`.claude/al_gamma_yield_check/`).

1. **`RESULT.md`** — the comparison, with a table of: thin-disk estimate,
   improved analytic estimate, Geant4 truth, and the ratio between them.
2. **A ranked breakdown of what the thin-disk approximation gets wrong**, each
   with a numerical factor and a sign. This is the actual question being asked.
   State plainly whether thin-disk is or is not adequate here.
3. **A clean simplified calculation** we can quote in a note — as few moving
   parts as possible while still being right to ~tens of percent.
4. If analytic and Geant4 still disagree after the corrections, **say so
   explicitly and name the leading suspects** rather than tuning the analytic
   model until it matches. A documented unresolved factor of 2 is more useful
   than a fudged agreement.

### Reporting conventions
- Rates in **per day** (× 1.929×10⁴ pulses/day) and **per pulse**, both.
- State the energy window on every number — thermal (<2 eV) and <1 keV are the
  two that get used, and the suspect calc used <1 eV. Do not mix them silently.
- Flag any nuclear-data value you could not verify against ENDF.

## 8. Related context

- `.claude/al_pair_background/VERDICT.md` — why this background matters:
  6×10⁵ MM-reconstructed pairs/day, ~10⁷:1 over X17 in the angular window,
  beatable only by a total-energy cut needing σ(p)/p ≲ 30 %.
- `docs/he3_self_shielding_note.md` — the He-3 optical-thickness argument.
- `CAMPAIGN_STATUS.md` — normalization conventions and current campaign state.
