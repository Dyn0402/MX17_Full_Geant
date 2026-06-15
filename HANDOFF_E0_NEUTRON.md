# Handoff: neutron-histories + E0 pair-channel work

**Last session:** 2026-06-15 · **Branch:** `main` (push before/after picking up).
For a **new Claude session**, possibly on another machine. Two collaborator
questions from Alberto drove this work; both have written-up answers, and the
E0 thread has a clear next step.

---

## 0. TL;DR — where things stand

- **Task 1 (neutron histories): DONE.** Answered + LaTeX note + event displays
  generated from a real lxplus Geant4 run. Nothing pending unless we want the
  optional E-at-capture scatter.
- **Task 2 (E0 pair channel): physics + final-metric DONE (2026-06-15).** The
  Formation/Decay note and the **final-metric note** (`e0_final_metric_note.tex`,
  `make_e0_final_metric.py`) are built. **Key result:** the ~10 µs flash readout
  only reaches E_n ≳ 20 keV, but E0 pairs are sub-keV → recorded E0 ≈ 0 in the
  flash window; even with a thermal trigger E0 gives ≤3 IPC bg pairs/day and
  ≤0.06 X17/day. **E0 is a model-completeness fix, not a game-changer.** Also
  settled: E0 *can* emit a massive vector/scalar X17 (not pseudoscalar/axial).
  See §4.
- The one open *physics* unknown is still the E0 strength `f = σ_E0/σ_M1`
  (bracketed `1e-3–1e-2`); pinning it needs a ⁴He R-matrix (a theory-contact
  ask, not a coding task). Remaining: discuss with Alberto (R-matrix; α_IPC
  reconciliation; whether a thermal trigger is worth scoping).

---

## 1. Environment (READ FIRST — cross-machine)

This work is **mostly local analysis + LaTeX**; only Task 1's sim ran on lxplus.

- **Python:** always `/home/dylan/PycharmProjects/nTof_x17/venv/bin/python`
  (has numpy, matplotlib, uproot, **h5py**). On a *different machine* this venv
  won't exist — recreate one with those four packages, or adjust the path.
- **LaTeX:** local `pdflatex` works (no `mdframed`/`amssymb`-beyond-base issues
  encountered; the notes only use base packages + `amssymb`). Build a note with
  `pdflatex <name>.tex` twice (for refs). `.aux` files are gitignored.
- **Geant4: NOT installed locally** — only ROOT 6.30 is. The sim runs on
  **lxplus**: `ssh lxplus` works from this machine (Kerberos + GSSAPI; if no
  ticket, `kinit dneff@CERN.CH`). Source tree on lxplus:
  `/afs/cern.ch/work/d/dneff/git/MX17_Full_Geant`. Build:
  `source scripts/setup_lxplus.sh && cd build && make -jN`. Deploy local C++
  changes by `scp` (verify `git diff <lxplus_HEAD> <local_HEAD> -- <files>` is
  empty so it's a clean overlay). **On another machine without the CERN
  Kerberos/SSH setup, lxplus access must be re-established first.**
- **BrIcc** is installed on lxplus at `/afs/cern.ch/work/d/dneff/tools/BrIcc`
  (`export BrIccHome=$PWD; ./bricc`). **Its tables stop at 6 MeV**, so it is
  useless at our 20.6 MeV — do not rely on it for the coefficient.

See memory `env_build_run.md` and `project_neutron_e0_questions.md` for more.

---

## 2. Task 1 — neutron histories (DONE)

**Question (Alberto):** are neutron scattering / slowing-down / capture in the
He-3 cell simulated, and could in-gas moderation enhance (or cause us to
undercount) the ⁴He* excitation rate?

**Answer:** Yes, fully (HP/G4NDL data, no neutron tracking cut). Counting
excitations via the (n,p)t channel, the full-transport rate per decade sits
at/below a single-pass estimate → not undercounting. Moderation is real
(75% E-loss per collision) but a tens-of-% perturbation in this thin 2 cm-bore
column (elastic MFP ≈ 24 cm vs 6 cm path).

**Deliverables (committed):**
- LaTeX: `docs/report/neutron_histories_note.tex` → `.pdf` (4 pp).
- Markdown: `docs/neutron_histories/neutron_histories_note.md`.
- Figures: `fig_excitation_fraction` (decade budget vs single-pass) and
  `fig_neutron_eventdisplay` (real Geant4 neutron paths at thermal/keV/MeV over
  the capsule cross-section, coloured by KE-fraction = slowing-down).
- **Sim feature added & built on lxplus:** `--trajdump [N]` flag dumps per-step
  neutron(+secondary) trajectories (`SteppingAction::DumpTrajectoryStep`,
  per-thread CSV). Run recipe: `docs/neutron_histories/event_display_RUN.md`.
  Plotter: `scripts/make_neutron_event_display.py`.

**Only optional remainder:** a higher-stats *E-at-capture vs E-incident*
scatter from the full campaign's terminal data on EOS (noted in the note).

---

## 3. Task 2 — E0 pair channel (IN PROGRESS)

**Question (Alberto):** the ⁴He 0⁺→0⁺ E0 transition makes e⁺e⁻ pairs by IPC
(γ-forbidden); is it in the sim, and how many IPC pairs per pulse vs energy?

### The physics, settled (all documented)
The problem splits into **Formation** (which ⁴He* state, set by angular
momentum) × **Decay** (how it reaches the 0⁺ g.s., set by selection rules),
linked by Jπ. Key facts:
- **E0 is missing from ENDF, Geant4, AND our generator** (it's γ-dark). Our
  sim's pairs are all hand-injected by `X17PrimaryGenerator`.
- **Two pair channels, opposite energy ends:** M1/E1 (0.35% IPC of a real γ,
  rides the rising σ_nγ → **peaks at MeV**, ~2200 pairs/day, computable now);
  E0 (100% pairs, s-wave + subthreshold-0⁺ boost → **piles up at sub-keV**).
- **0⁺→0⁺ goes ONLY by E0** (M1 forbidden too; same 0→0 rule). The M1 is a
  *different* capture channel (³S₁→1⁺).
- **No low-lying 1⁺** (lowest ~28 MeV) → M1 capture is **direct/continuum**, not
  resonant. In fact all ⁴He levels are broad (Γ/|E_r|~1.2–1.9) → capture is
  **mostly direct**; only the subthreshold 0⁺ gives a (broad) ~10× low-E boost
  to the E0.
- Only **three** states make g.s. pairs: 0⁺ (E0), 1⁺ (M1), 1⁻ (E1). 0⁻ and 2⁻
  are spectators (0⁻→0⁺ forbidden for everything; 2⁻→0⁺ is M2, too slow).

### Deliverables (committed)
- **LaTeX (the main note):** `docs/report/e0_pair_channel_note.tex` → `.pdf`
  (6 pp), structured Part 0 summary / Part 1 Formation / Part 2 Decay / combine.
- **Markdown spine:** `docs/e0_branch/formation_and_decay.md`. Plus
  `capture_breakdown.md`, `doorway_states_note.md`, `e0_transition_explainer.md`,
  `ipc_estimation_method.md` (has the literature dig + **References**),
  `ipc_roadmap.md`.

**Part 1 expanded (2026-06-15)** — `doorway_states_note.md` rewritten as the full
Part-1 reference (TUNL data provenance; continuum-vs-resonance physics; the
formation decomposition `S=(1/v)·g_J·P_ℓ·M`; the relative-probability recipe +
R-matrix/AZURE2 pointer; isospin-forbidden E1 → why M1 dominates low-E radiative
capture, with Wolfs/Wervelman refs). **Two figure corrections:**
(1) `fig_he4_levels` now has a y-axis break 4→16 MeV; (2) `fig_doorway_formation`
**rewritten** — the old "1⁺ dominates by ~10×" was an ad-hoc-normalisation
artifact (flat `F=4/Γ_ref` for the direct 1⁺). Honest result: 1⁺ and 0⁺ are
COMPARABLE at thermal (`S(1⁺)/S(0⁺)=3/boost_0⁺ ≈ 0.5–1`; the ×3 spin weight is
cancelled by the 0⁺ sub-threshold boost), gap opens only toward MeV. Report-note
figure copies (`docs/report/figs/fig_e0_{levels,doorway,directres}.*`) refreshed.
NB: `e0_pair_channel_note.tex` Part-1 *text* not yet updated for the isospin
point / corrected comparison (figures are current).

**Cross-section note (2026-06-15)** — new standalone note
`docs/report/e0_cross_section_note.tex` → `.pdf` (3 pp), the clearest framing so
far: σ = σ_form × (Γ_X/Γ_tot), then (1) **capture cross sections**
(`make_e0_cross_section_fig.py` → `fig_e0_xsec`): measured σ_nγ (M1+E1, ENDF,
55 µb thermal) vs estimated E0 band (f=1e-3–1e-2), E0 sitting 4–5 decades below;
(2) **IPC yield = σ × pair fraction** (`make_e0_ipc_yield_fig.py` →
`fig_e0_ipcyield`): ×α_IPC≈3.5e-3 for M1+E1, ×1 for E0 → the gap collapses, pair
yields **comparable at thermal** (E0/M1E1 = f/α_IPC ≈ 0.3–3), E0 sub-keV / M1E1
MeV. Both read `data/He3.h5`. Uncertainty table: f dominates (~1 decade); σ_nγ
firm, α_IPC ~×1.8. Dylan's reaction: clearest plot yet.
Now 5 pp — extended with **(3) X17 production, two scenarios**
(`make_e0_x17_production_fig.py` → `fig_e0_x17`): X17 = IPC yield × BR_X17(0.025);
scenario 1 (vector/scalar X17) = M1+E1 + E0, scenario 2 (pseudoscalar/axial) =
M1+E1 only (the E0/monopole selection rule is the discriminator). **Critical
window stats (X17 produced/day):** 0.2–2 MeV high-E window = **37.4/day**
(22.5/day at table α; E0 adds <0.004 → identical in both scenarios); full sub-keV
window = **0.15/day** (scenario 2) → **0.18–0.42/day** (scenario 1, f band). Key
conclusion: the X17 rate lives in the high-E window and is X17-type-independent;
the sub-keV is starved either way (≤0.4/day even for vector/scalar X17, ~100×
below the high-E window) — E0 adds a *distinct sub-keV signature*, not
*statistics*. Recorded ≈ ×0.196; sub-keV needs a thermal trigger.
Now 6 pp — added **(4) "From cross section to statistics" bridge**
(`make_e0_luminosity_fig.py` → `fig_e0_lumi`, 4 panels): N_X/day = φ × f_breakup
× (σ_X/σ_np) × ppd. Key clarification (Dylan asked): for E0 (not in Geant4),
N_E0 = N_(n,p) × σ_E0/σ_np — the (n,p)t breakup is the luminosity monitor, and
**opacity cancels in the ratio** (the (n,p) absorption attenuates both channels
equally, so N_X = N_np σ_X/σ_np holds exactly even sub-keV). M1+E1 uses the
direct Geant4 (n,γ) count. Panels: (a) EAR2 flux, (b) breakup fraction
opaque(70%)→thin(4%), (c) effective luminosity L=N_np/σ_np, (d) breakups per X17
= σ_np/σ_X17, falling 1e12 (sub-keV) → 1e8 (MeV) — the quantitative reason the
MeV window is favourable.

**Folded capstone (2026-06-15)** — note now 6 pp, closes with
`make_e0_folded_fig.py` → `fig_e0_folded`: 2 panels (X17 left, IPC right),
per-day-per-decade, PRODUCED and RECORDED (×MM acceptance 0.196 X17 / 0.236 IPC),
scenario-1 band (+E0) vs scenario-2 (M1+E1). On-plot text boxes give window sums.
**Window totals (produced → recorded):** X17 0.2–2 MeV = 37.4→7.3/day;
X17 sub-keV = 0.15→0.03 (sc2) / 0.18–0.42→0.035–0.083 (sc1). IPC 0.2–2 MeV =
1497→353/day; IPC sub-keV = 6.1→1.4 (sc2) / 7.2–16.8→1.7–4.0 (sc1). The
`e0_cross_section_note.tex` is now the cleanest end-to-end note: σ=form×branch →
capture xsec → IPC yield → X17 (2 scenarios) → flux/luminosity bridge → folded
per-day capstone → uncertainties. Five figure scripts: make_e0_{cross_section,
ipc_yield,x17_production,luminosity,folded}_fig.py. NOT yet committed to git.
- **Figures + scripts** (all in `scripts/`, output to `docs/e0_branch/figs/`):
  | script | figure | shows |
  |---|---|---|
  | `make_he4_levels_fig.py` | `fig_he4_levels` | ⁴He level diagram, TUNL A=4 |
  | `calc_doorway_formation.py` | `fig_doorway_formation` | s-wave (0⁺,1⁺) vs p-wave (0⁻,2⁻) formation vs E_n |
  | `calc_direct_vs_resonant.py` | `fig_direct_vs_resonant` | which doorways are direct vs resonant |
  | `make_capture_breakdown_fig.py` | `fig_capture_breakdown` | E0 vs M1 channel schematic |
  | `make_ipc_vs_energy_fig.py` | `fig_ipc_vs_energy` | term-1 (M1/E1) pairs/day vs E_n |
  | `make_e0_pair_yield_fig.py` | `fig_e0_pair_yield` | E0 vs M1/E1, opposite energy ends |
  | `plot_he3_xs.py` | `fig_he3_xs` | ENDF (n,p)/(n,γ)/ratio (needs h5py) |

### The E0 strength `f` (the open unknown)
`f = σ_E0/σ_M1` is bracketed **1e-3–1e-2** by a transparent argument (statistical
¼ ¹S₀ weight floor; E0 has no real-photon mode ceiling) → E0 pairs ≈ 0.3–3× the
M1 pairs at low E_n. No literature value exists; the ⁴He(e,e′) monopole data fix
the *decay-side* matrix element but the capture rate needs a **⁴He multi-level
R-matrix** (Hale's evaluated system) normalised by that strength. This is a
theory-contact ask, the rate-limiter for a precise number. Details +
all refs: `docs/e0_branch/ipc_estimation_method.md`.

---

## 4. NEXT STEP — final-metric note (DONE 2026-06-15)

**The "final-metric" note is built.** `N_rec(E_n) = [N_M1/E1 + N_E0(f)] ×
acceptance × TOF-window`, summed over decades × pulses/day. Deliverables:
- **Script:** `scripts/make_e0_final_metric.py` → `analysis/e0/final_metric.json`
  + `docs/report/figs/fig_e0_final_metric.{pdf,png}` (two panels: IPC background
  and X17 signal, produced vs recorded, with the flash/thermal trigger regions).
- **LaTeX:** `docs/report/e0_final_metric_note.tex` → `.pdf` (4 pp).
- Term 1 = direct G4 (n,γ) counts × α_IPC (reproduces the committed mev_note:
  32.7 X17/day produced, 6.3/day recorded ≈190/30-day at table α=2.1e-3). Term 2
  = the existing sub-threshold-0⁺ E0 construction over f=1e-3–1e-2. Acceptance =
  flat MM-double 19.6% X17 / 23.6% IPC (angular note).

**The decisive finding — the trigger window kills recorded E0.** At L=19.5 m the
~10 µs flash readout reaches only **E_n ≳ 20 keV** (TOF=10 µs ⇔ 19.9 keV); the
E0 pairs are **sub-keV** (TOF >45 µs) → arrive after the window. So:
- **Flash readout:** ~500 IPC bg pairs/day + ~6–10 X17/day recorded, **all
  M1/E1**; recorded E0 ≈ 0. Baseline run is unaffected by E0.
- **+ thermal trigger** (the only way to reach E0): even at f=1e-2, E0 gives
  only ≤3 IPC bg pairs/day and ≤0.06 X17/day — because the E0 capture inherits
  the σ_nγ/σ_np ~1e-8 sub-keV suppression × f. S/B (vs ~3e5 (n,p)t/pulse), not
  rate, is the obstacle. **E0 is a completeness fix, not a game-changer.**

**E0 → X17 settled (Dylan's question):** a massive X17 is NOT γ-dark to a 0⁺→0⁺
monopole (the photon's masslessness is what forbids it). **Vector(1⁻)/scalar(0⁺)
X17 CAN be emitted in E0; pseudoscalar(0⁻)/axial(1⁺) cannot.** The E0→X17
branching is a separate unknown (not the M1's 2.5%); benchmarked at 2.5% in the
table. Rate still tiny because the E0 capture is tiny.

**Trigger framing (Dylan, 2026-06-15):** now considering BOTH the 10 µs flash
readout AND a separate thermal-neutron trigger — the latter is what would give
sub-keV statistics where E0 (and any E0→X17) lives.

**Still parallel / when possible:** raise the E0 R-matrix question with Alberto
or a ⁴He-theory contact to replace the `f` bracket with a number. Also reconcile
α_IPC = 3.5e-3 (our anchor) vs 2.1e-3 (Alberto's table) with him.

---

## 5. Pointers
- Memory (persists across sessions, this machine): `MEMORY.md` index →
  `project_neutron_e0_questions.md` (the full blow-by-blow),
  `env_build_run.md`, `project_mx17_full_sim.md`. On another machine the memory
  dir won't be present — this handoff is the portable substitute.
- All work is on `main`, committed. 10 commits this session
  (`08c9d88`…`30b1d85`). Older campaign context: `CAMPAIGN_STATUS.md`,
  `HANDOFF_FULL_SIM.md`, `docs/report/{thermal,mev,angular}_note.tex`.
