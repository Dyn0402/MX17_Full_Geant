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
- **Task 2 (E0 pair channel): IN PROGRESS, well advanced.** Physics fully worked
  out and documented (LaTeX note rebuilt on a Formation/Decay spine). **The next
  concrete step is the "final-metric" calculation:** combine the formation
  shapes × decay pair-fractions × detector acceptance → **X17 / IPC pairs per
  pulse (per day) vs neutron energy, including the E0 channel.**
- The one open *physics* unknown is the E0 strength `f = σ_E0/σ_M1` (bracketed
  `1e-3–1e-2`); pinning it needs a ⁴He R-matrix (a theory-contact ask, not a
  coding task).

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

## 4. NEXT STEP (where to start)

**Build the "final-metric" note:** combine
`N_pairs(E_n) = Σ_Jπ [formation shape] × [pair fraction]`, fold in the absolute
normalisations (campaign (n,p)t rate, α_IPC≈3.5e-3, and the bracketed f for E0)
and the detector acceptance, to produce **X17 / IPC pairs per pulse and per day
vs neutron energy, with the E0 channel included** (and shown as a band over the
f range). Inputs already on disk:
- campaign rates: `analysis/mev/mev_rates.json` (per-decade (n,p)t and direct
  (n,γ) counts, `pulses_per_day`, weights).
- term-1 machinery: `scripts/make_ipc_vs_energy_fig.py`.
- E0 shape machinery: `scripts/make_e0_pair_yield_fig.py`.
- acceptance: pairs-acceptance results live in `analysis/pairs_v2/` and the
  angular-resolution work (`docs/angular_resolution/`); check how X17/day was
  computed in the existing MeV note (`docs/report/mev_note.tex`) to stay
  consistent.

**Parallel / when possible:** raise the E0 R-matrix question with Alberto or a
⁴He-theory contact to replace the `f` bracket with a number. Also reconcile
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
