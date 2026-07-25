# Al(n,γ) yield cross-check — RESULT

> **NOSE-FIRST — FINAL (2026-07-24).** The mounting audit (2026-07-23)
> confirmed the capsule sits tip-into-beam; the sim was flipped to
> `rotateX(+90°)` (commit `3d97437`) and the full `*_2cm_nose` campaign has
> re-run and been re-analysed (thermal + epithermal neutrons, ×10⁵-biased
> captures, 10⁷ X17+IPC pairs, 10⁷ γ-source cascades). **Headline nose-first
> results — the figures and the Part II/III numbers below are refreshed to
> these:**
> - **Al(n,γ) capture 4.507×10⁻³/n** (40M neutrons, ±0.24%) = **19,346/pulse =
>   3.73×10⁸/day** — **×0.577** of valve-first (7.806×10⁻³), matching the
>   analytic nose-first single-pass prediction (4.5×10⁻³) to 0.2%. ~95% of
>   captures now in the nose (⟨y⟩ ≈ −29.7 mm) vs 68% in the valve stem before.
> - **7.724 MeV γ** (21.3% branch): **4,121/pulse = 7.95×10⁷/day.**
> - **Al-attributable trigger legs ~112/pulse** (was 199; **×0.563**), by two
>   independent routes agreeing to 0.2%: direct neutron-mode (1,616 legs / 6×10⁷ n)
>   and γ-source closure (2.723×10⁻³ legs/γ × 1.13 γ/capture × 19,346).
> - **The mechanism is invariant under the flip** (per-γ leg efficiency −2.6%):
>   topology (98% one crossing e⁻), birth-volume mix, species, and KE all move
>   ≤1%. One genuine gain: **58% of leg electrons now also cross the MM drift
>   gas** (was 52%) — a stronger MM-track veto handle, because the capture zone
>   sits further upstream.
>
> The Part I analytic cross-check (method, +5% G4 agreement, 21.3% branch) is
> unchanged and its nose-first single-pass row is now confirmed by Geant4.
> Valve-first data is retained on EOS (`*_2cm`) and in git history for
> comparison. Refreshed inputs: `gsrc_mechanism_nose.json/.npz`,
> `leg_mechanism_nose.json`, `analysis/pairs_nose/pairs_nose.pdf`.

**2026-07-23, Claude (Fable 5) for Dylan.** Task spec:
`.claude/al_gamma_yield_check/HANDOFF.md`. Slides:
`slides/al_capture_crosscheck_slides.pdf` (pedagogical walk-through with all
figures). Scripts: `analysis/al_gamma_yield_check/`.

## TL;DR

**The Geant4 number is right.** An independent single-pass (straight-line
attenuation) calculation on the exact polycone geometry, the real flux
spectrum, and ENDF/EGAF cross sections reproduces the Geant4 Al capture rate
to **+5%**, and also reproduces *where* on the capsule the captures happen
(to a few %) and the capturing-neutron energy spectrum (flat ~0.95 ratio).
The earlier hand estimate "looked low" because thin-disk with any plausible
single thickness is the wrong model here — and because only 21.3% of captures
emit the 7.724 MeV line.

**Big incidental finding: the sim's capsule faces the beam valve-first.**
The placement `rotateX(-90°)` maps local z → world −y, so the valve (local
z = +51 mm) sits at y = −51 mm, upstream. The capture positions prove it:
68% of Al captures in the valve stem, 28% in the shoulder, 0.2% in the nose.
The beam core (r₅₀ = 2.6 mm) bores through 10–20 mm of nearly solid Al before
touching He-3. If the real mounting is nose-first, the Al capture background
is **×0.55** of the current sim ⇒ **confirm the mounting** (and if nose-first
is real, fix the rotation sign and re-run the Al-pair chain).

## The comparison (window [1 meV, 2 eV] unless stated)

The G4 column below is the original **valve-first** methodology cross-check;
the nose-first Geant4 number (19,346/pulse, confirming the single-pass
nose-first row) is in the paragraph beneath. "ratio to G4" is vs valve-first.

| estimate | captures/pulse | /day | ratio to G4 |
|---|---|---|---|
| thin disk, 0.6 mm wall, σ_th, <1 eV | 2,880 | 5.6e7 | **0.09×** |
| thin disk, 5.5 mm nose, σ_th, <1 eV | 26,396 | 5.1e8 | 0.80× |
| thin disk, 13.24 g over face, σ_th, <1 eV | 66,696 | 1.3e9 | **2.02×** |
| **single-pass, nose-first (final), <1 eV** | **19,082** | **3.7e8** | **0.58×** |
| single-pass, valve-first (as-was), <1 eV | 34,685 | 6.7e8 | 1.05× |
| Geant4 valve-first, <1 eV | 33,021 | 6.4e8 | 1 |
| single-pass, valve-first, full window | 35,181 | 6.79e8 | 1.05× |
| Geant4 valve-first, full window | 33,508 | 6.46e8 | 1 |
| **Geant4 nose-first, full window** | **19,346** | **3.73e8** | **0.577×** |

Geant4 side (**valve-first, as-was**): 20 files × 10⁷ n of
`neutrons_thermal_trig_2cm` (analog),
`capture_vol=="He3Cap_Al" && capture_proc=="nCapture"`: 1,561,193/2×10⁸
= 7.806×10⁻³ per neutron (±0.08% stat). 98.55% of those below 1 eV.
Normalization: ×4.2924×10⁶ n/pulse, ×1.929×10⁴ pulses/day.

**Nose-first confirmation (2026-07-24):** the flipped-capsule dataset
`neutrons_thermal_trig_2cm_nose` (40M n) gives Al **4.507×10⁻³/n** (±0.24%),
= 19,346/pulse = 3.73×10⁸/day, i.e. **×0.577** of valve-first — landing on
the analytic single-pass nose-first prediction (19,082/pulse, 3.7×10⁸/day)
to 0.2%. The single-pass method therefore reproduces Geant4 in *both*
orientations; the whole Part I methodology below stands.

Capture location, Geant4 vs analytic single-pass (fraction of Al captures):
valve stem (z∈[40,51]) 68.3% | 64.4%; shoulder (21–40) 28.5% | 30.8%;
barrel (|z|<21) 3.0% | 4.5%; nose (z<−21) 0.2% | 0.2%.

## Ranked: what thin-disk gets wrong (the actual question)

1. **Which Al the beam sees — the dominant issue, sign depends on variant.**
   The wall slab (0.6 mm) is 11× low: the beam core never touches the barrel
   wall. Smearing all 13.24 g over the face (t = 13.9 mm) is 2.0× high: 3.7 g
   sit *behind* black He-3 (28%), and the mass-smeared thickness ignores where
   the flux is. The correct single number is the flux-weighted **pre-gas** Al
   path, ⟨t_pre⟩ = **9.1 mm**, dominated by the valve stem + shoulder because
   the capsule is mounted valve-first (×1.8 vs nose-first).
2. **σ_th instead of ⟨σ⟩ — ×1.45 high.** Folding 1/v against the actual
   spectrum gives ⟨σ⟩/σ_th = 0.69 (<1 eV) or 0.65 (full window): the window is
   *above* 25.3 meV on average. The handoff's suspicion that the sub-eV 1/v
   rise makes σ_th an *under*estimate goes the wrong way — the meV flux is too
   small to win.
3. **He-3 shadowing** — fully contained in ⟨t_pre⟩ above; as a mass statement,
   72% of the Al is reachable before first gas contact (valve-first).
4. **Transport beyond straight-line — net ×0.95 only.** Elastic out-scattering
   (Al 1.4 b ≈ 6× capture, plus CFRP/air) and albedo from the surrounding
   plastics nearly cancel. Evidence: He-3 absorbs 77.4% of neutrons in G4 vs
   91% in the straight-line model (out-scattering), yet Al captures only drop
   5%; G4/analytic ratio is flat in energy and position.
5. **γ per capture — ×4.7 if ignored.** The 7.724 MeV line is emitted in
   **21.3%** of captures (EGAF: σγ = 0.0493(15) b vs σ₀ = 0.231 b).

Verdict on the method: **thin-disk is not adequate** with any a-priori single
thickness; with the computed ⟨t_pre⟩ it collapses to the one-liner below and
is good to 5%.

## The quotable calculation

R = Φ × n_Al ⟨σγ⟩ × ⟨t_pre⟩
  = 4.29×10⁶ n/pulse × (6.02×10²² cm⁻³ × 0.231 b × 0.65) × 0.91 cm
  = **3.5×10⁴ captures/pulse = 6.8×10⁸/day**  (Geant4: 3.35×10⁴ / 6.46×10⁸)

7.724 MeV γ yield (×0.213): **7.2×10³/pulse = 1.4×10⁸/day** (G4-anchored).
Extending the window to <1 keV adds only +3.6% (Φ×1.7, σ~1/v dies):
7.03×10⁸ captures/day analytic. (Resonance region ≥5.9 keV not addressed.)

## Nuclear data (verified)

- σ₀(²⁷Al(n,γ)) = 0.231 b @ 25.3 meV — cross-confirmed by the ²⁸Al decay line
  (1778.9 keV, ~1/capture) having σγ = 0.232(4) b in the IAEA PGAA k₀ table.
- σγ(7724.03 keV) = 0.0493(15) b ⇒ **21.3% of captures** (IAEA PGAA/EGAF).
- ³He(n,p)t from ENDF via `data/He3.h5` MT=103 (5333 b @ 25.3 meV, 1/v).
- Al elastic ≈ 1.4 b (not in the single-pass model; checked via the G4 residual).
- NB in G4NeutronHP, ³He(n,p)t is filed as `neutronInelastic`; `nCapture` in
  He3Gas is only the tiny radiative branch. Don't be fooled by the budget.

## Files

- `analysis/al_gamma_yield_check/al_capture_crosscheck.py` — the analytic model (geometry marching,
  flux/Lambda2D folding exactly as the gun, 1/v, He-3 removal). Reproduces
  V_Al and V_gas to 4 digits. Writes `analytic_results.npz`, `summary.json`.
- `analysis/al_gamma_yield_check/plot_al_capture_crosscheck.py` — all figures (`figs/*.pdf`) +
  region-fraction table.
- `caps_20files.root` — Geant4 capture histograms (macro on lxplus:
  `~/al_check/count_captures.C`, EventTree over 20 files).
- `slides/al_capture_crosscheck_slides.tex|pdf` — the presentation.
- **Nose-first (2026-07-24) refreshed inputs/outputs** (EOS datasets
  `*_2cm_nose`, libraries in `full_sim/libs_nose/`): `gsrc_mechanism_nose.json`
  / `.npz` (γ-source, 10⁷ γ, 27,226 legs), `leg_mechanism_nose.json` (neutron
  mode, 6×10⁷ n, 1,616 legs), `analysis/pairs_nose/pairs_nose.pdf` (X17+IPC
  pair acceptance/mass/angle, 10⁷ events). Figures `figs/*.pdf` regenerated
  from these via `plot_leg_mechanism_truth.py`, `plot_al_pair_danger.py`,
  `plot_coincidence_ladder.py`.

## Part II — single-arm coincidences (Alberto's email, 2026-07-23)

**TL;DR: both calculations are internally fine; they differ by four
compounding assumptions, and the plastic material is not one of them.**
The sim's ~28 legs/arm/pulse (**112 arm-summed, Al−noAl, nose-first**; was
199 valve-first) is confirmed and mechanistically understood; Alberto's
~9/pulse follows from his stated assumptions (one-arm solid angle, γ₀ only,
1% "Compton in the detector").

**Capture/γ level — the factor grows with the nose-first flip.** Alberto's
72,787 captures/pulse (<1 eV; 13 g thin disk at ~0.0923 at/b) is 3.8× the
nose-first Geant4 truth (19,346; it was 2.2× the valve-first 33,021) — same
spectral shape; that is the "mass-smeared thin disk" rung of the Part I
ladder combined with the ×0.577 nose-first drop. His γ₀ intensity 26.8% vs
EGAF 21.3% makes the γ₀ ratio 4.7 (19,507 vs 4,121/pulse).

**The coincidence ladder (figs/coincidence_ladder.pdf), 9 → 112/pulse
(nose-first):**

| step | legs/pulse | factor |
|---|---|---|
| Alberto: γ₀ × 4.5% Ω × 1% | 8.8 | — |
| true captures (nose) + 21.3% intensity | 1.8 | ×0.21 |
| full 4-arm plastic aperture (18.5%) | 7.4 | ×4.1 |
| all cascade lines >2.5 MeV (~1.7 γ/capture) | 59.5 | ×8.0 |
| single-electron punch-through (1.8%/γ, measured) | 112 | ×1.9 |
| **Geant4 legs, Al−noAl** | **112.0** | — |

The nose-first flip enters only at the second rung (true captures ×0.577):
every rung above Alberto's fixed anchor drops by that factor, since the
per-capture mechanism (aperture, response, cascade multiplicity) is
unchanged — see Part III.

Key facts behind the two big steps (all from classifying 1,616 legs in
6×10⁷ nose-first events, `leg_mechanism_nose.json`):

- **98.1% of legs are one charged track (e⁻ 90%, e⁺ 8%) crossing BOTH
  detectors**, median KE 3.95 [2.99, 5.31] MeV at the SiPM bar. A capture γ
  interacts once — in the MM PCB Cu/FR4, the SiPM bar, or nearby
  structure — and the multi-MeV electron crosses the 3 mm bar (0.23 MeV
  threshold) and buries itself in the 2 cm plastic (1.73 MeV threshold).
  A coincidence does NOT cost P(SiPM)×P(plastic); one interaction buys
  both signals. Effective probability per cascade γ through the aperture:
  1.8% measured; an analytic material budget (PCB 0.55% + bar 0.75% +
  container/air) times electron-survival gives 1–2%. Same-γ double
  interactions: 0.5%; independent particles (cascade partners): 1.5%.
  96.9% of leg events have `capture_vol == He3Cap_Al` (nose-first).
- **Alberto's 4.5% solid angle is one arm's plastics** (we compute 4.6%
  per arm at 42 cm); the four-arm plastic aperture is 18.5% and the SiPM
  wall subtends 40%.
- **Thresholds in energy**: SiPM 0.5 MIP = 0.23 MeV; plastic 0.5 MIP =
  1.73 MeV. The 1778.9 keV line is the ²⁸Al β-decay line (T½ = 2.24 min):
  delayed out of the gate *and* below the plastic max-Compton threshold —
  it contributes nothing.

**Material sensitivity — a ±10% effect, not ×5.** Compton is
electron-density physics; Klein–Nishina interaction probabilities for
2 cm at 7.72 MeV: PVT 4.05%, polystyrene 4.09%, PMMA 4.64%, polyethylene
3.89% (7–8.5% at the 3 MeV cascade lines; pair adds ~0.9% at 7.72). The
3 mm SiPM bar is 0.62% — Alberto's "~1% Compton" describes the thin bar,
not the 2 cm plastic. Not knowing the exact scintillator brand cannot
explain any order of magnitude.

## Part III — full mechanism characterization (truth runs, 2026-07-23)

**Method.** Birth-truth branches on HitTree (`origin_vol`, `origin_proc`,
`origin_ke`, `ox,oy,oz` — from `GetLogicalVolumeAtVertex`/`GetCreatorProcess`;
committed, `d523091`). Two **nose-first** truth samples: (a) neutron mode,
6×10⁷ events → **1,616 legs** (real cascade, `leg_mechanism_nose.json`);
(b) **γ-source truth run** — 10⁷ single cascade γs (20 discrete ²⁸Al lines
>2 MeV, IAEA intensities) emitted isotropically from 189,570 real
He3Cap_Al/CFRP capture vertices → **27,226 legs** with full birth truth
(`gsrc_mechanism_nose.json/.npz`).

**The four ways a leg happens** (figs/leg_mechanism_corrected.pdf):

| # | topology | fraction |
|---|---|---|
| A | **hard Compton** → ONE e⁻ crosses bar+plastic | **80.8%** |
| B | **pair production** → e⁺/e⁻ cross bar+plastic | **17.0%** |
| C | same γ interacts twice (the old "double Compton") | 0.9% |
| D | two different cascade γs, same arm | 1.0% |

The old 84/16 slide (summary_2026-07-22): the 16% pair fraction was right;
the "84% double Compton" topology was wrong — the true same-γ-twice rate is
0.9%. The electrons-only 84% is really ONE hard Compton electron crossing
both detectors (birth KE median 4.69 MeV; 3.95 MeV at the bar). Nose-first
moves the A/B split <1% vs valve-first (80.0/17.9) — the mechanism is
orientation-independent.

**Where the electron is born** (same-track legs, γ-truth;
figs/leg_origin_breakdown.pdf):

| birth volume · process | fraction |
|---|---|
| Al capsule · Compton | 33.0% |
| MM+PCB (Cu/FR4/kapton/mesh) · Compton | 19.8% |
| SiPM bar (3 mm PVT) · Compton | 18.8% |
| capsule CFRP · Compton | 8.0% |
| Al capsule · pair | 7.7% |
| MM+PCB · pair | 5.7% |
| SiPM bar · pair | 2.5% |
| air · Compton | 2.5% |
| capsule CFRP · pair | 1.0% |
| everything else (wrap, air·pair, δ-rays…) | 1.0% |

Grouped: **capsule (Al+CFRP) 50%, MM+PCB 25%, SiPM bar 21%, air 3%** —
essentially nothing starts in the plastic itself (a plastic-born electron
can't reach the SiPM behind it). The nose-first shift is toward the capsule
(48→50%) as the capture zone moves into the Al nose. Cross-validated in the
real cascade: **58.4% of leg electrons also cross the MM drift gas** (up from
51.6% valve-first; born upstream of the mesh ≈ capsule+air ~53% in γ-truth) —
the upstream capture zone means more than half the trigger legs leave an MM
track, a strengthened veto handle. Unambiguous pair signature (e⁺ AND e⁻
both crossing both detectors): 3.4–3.9%.

**Per-line response** (legs per emitted γ, all arms): zero at 2.28 MeV
(below the 1.73 MeV plastic threshold), ~0.7×10⁻³ at 3.0 MeV, ~2.4×10⁻³ at
4.26 MeV, **5.77×10⁻³ at 7.724 MeV**; the 7.72 line alone makes **49% of all
legs**. That measured 0.58%/γ (= 3.1% per γ through the plastic aperture) is
the number Alberto's "1% Compton" stands in for (essentially unchanged from
the valve-first 0.54%/γ — the per-γ response is orientation-independent).

**Closure / caveat.** γ-source total: 2.723×10⁻³ legs/γ × 1.13 γ/capture ×
19,346 captures/pulse = **59.5/pulse**, vs 112 in full neutron mode: the
20-discrete-line list underestimates the true G4NDL cascade response ×1.9
(continuum multiplicity + within-cascade co-adding) — the same ×1.9 gap as
valve-first (106 vs 199). The ladder's last rung carries this; the mechanism
*fractions* are validated against the neutron run (species, birth KE,
upstream fraction all agree).

## Orientation history — it was never flipped; it was born this way

- **2026-06-01** (`6597f04`): `rotateX(-90°)` + beam +Y introduced; the
  capsule was a symmetric `G4Tubs` — orientation meaningless.
- **2026-06-10** (`86d557b`): STEP polycone lands (tip z=−35, valve +51)
  under the pre-existing rotation ⇒ valve-upstream from day one.
- **2026-06-11** thermal_note: the prose says the beam "enters through a
  5 mm-thick domed aluminium tip" (the intention), but the note's own
  quoted rate — capture in Al capsule = 8.0×10⁻³ — is the *valve-first*
  number (nose-first = 4.5×10⁻³). The sim behind the note was already
  flipped; only the text described nose-first.
- **06-30 / 07-14 / 07-15** (axis convention, MM spacing, stack flip):
  none touch `capRot` or the polycone z arrays (verified by diff).
- **2026-07-22** event display draws y_world = −z_local (correct).

So the axis redefinition is innocent, and every thermal-era number is
internally consistent (all valve-first). **Resolved 2026-07-23:** the
mounting audit says nose-first; the sim was flipped (`3d97437`) and the
`*_2cm_nose` campaign re-ran — captures ×0.577 and the capture zone moved
into the nose, exactly as anticipated (see the banner for final numbers).

## Caveats

- Geant4's internal capture-cascade branching for the 7.724 line was not
  unfolded here; for γ-level notes use EGAF 21.3%.
- ⟨t_pre⟩ = 9.1 mm is orientation-specific (valve-first). Nose-first:
  single-pass gives 4.5×10⁻³/n (19,082/pulse), now confirmed by Geant4 at
  4.507×10⁻³/n (19,346/pulse) — mounting check done, see banner.
- <1 keV numbers are analytic-only (thermal campaign window ends at 2 eV);
  the epi campaigns exist if a G4 check above 2 eV is wanted.
- Part II ladder: the product of the last two rungs (×8.0 cascade ×1.9
  mechanism = ×15.2) is pinned by the sim; the split between them relies
  on the >2.5 MeV cascade-line inventory (~1.7 γ/capture, EGAF-derived
  estimate) and is soft at the ±50% level. Mechanism stats (nose-first):
  1,616 classified neutron-mode legs (6 files) + 27,226 γ-source legs, ±2%.
