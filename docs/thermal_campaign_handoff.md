# Handoff: neutron campaign + thermal-statistics analysis (2026-06-11)

One-page-ish summary of everything done on 2026-06-10/11, for a future
session picking this up cold. Companion docs:
[../CAMPAIGN_STATUS.md](../CAMPAIGN_STATUS.md) (live status),
[../PLAN_NEUTRON_CAMPAIGN.md](../PLAN_NEUTRON_CAMPAIGN.md) (original plan),
[he3_self_shielding_note.md](he3_self_shielding_note.md) (analytic argument),
[report/thermal_note.pdf](report/thermal_note.pdf) (the deliverable).

---

## The physics result (the thing that matters)

**The sub-keV X17 measurement is not feasible.** The collaboration rate
table (`/home/dylan/x17/calculation_tables/results_3He`) computed radiative
captures thin-target (column × σ_nγ), but below 1 keV the ³He gas is opaque
to (n,p)t — every beam neutron is absorbed regardless, so radiative capture
per beam neutron is capped at σ_nγ/σ_np ≈ 1.0×10⁻⁸ (thermal), rising to
~2×10⁻⁸ at 1 keV (ENDF/B-VIII.0; see `scripts/plot_he3_xs.py`).

Run B (10⁹ fully-transported Geant4 neutrons, EAR2 flux + Lambda2D
footprint) confirms this empirically. Integrated over E_n = 1 meV–1 keV,
with 1.929×10⁴ pulses/day (7×10¹² ppp):

| quantity | table /pulse | sim /pulse | table /day | sim /day |
|---|---|---|---|---|
| radiative captures | 5.8 | (0.5–1.1)×10⁻¹ | 1.1×10⁵ | (1.0–2.1)×10³ |
| IPC pairs | 1.21×10⁻² | (1.1–2.3)×10⁻⁴ | 233 | 2.2–4.4 |
| X17 (2.5% **of IPC**) | — | (2.9–5.7)×10⁻⁶ | (5.8) | **0.06–0.11** |

→ **~1 X17 produced every 9–18 days, before acceptance.** Factor 50–100
below the table (×53±14 from 15 direct (n,γ) events, ×106 from the
thermal-ratio floor). Sub-keV-anchored sensitivity: 3.0σ → 0.3–0.4σ.

**Likely resolution: move the ROI to 0.1–10 MeV** where the gas is thin,
the table rows are valid, and ~98% of the physical rate lives (the
commented-out 0.2–2 MeV row in `results_3He`: 5.2×10⁻² IPC/pulse). This
needs (a) E_n-dependent generator kinematics (E* ≈ 20.58 MeV + ¾E_n + CM
boost), (b) gamma-flash/TOF assessment at µs timescales. **Run B-full
(1 meV–100 MeV, 5×10⁸ events) is already on disk for the capture budget.**

### Pitfalls a future session must know

- **2.5% X17 branching is per IPC pair, not per radiative capture** (a
  factor-500 bug fixed on 06-11; commit 8d3453e).
- **³He(n,p)t appears as `neutronInelastic`, not `nCapture`** in Geant4 HP;
  `nCapture` in He3Gas = radiative capture = the X17 channel.
- The 2.1×10⁻³ IPC/capture coefficient and the 2.5% are *assumptions* —
  provenance to be confirmed with Alberto (multipolarity? E0?).
- Direct (n,γ) statistics: 15 events in 10⁹ — fine for the ×50–100
  conclusion, not a precision number.

### Secondary findings (backgrounds, from the same scan)

Per beam neutron < 1 keV: Al capsule nCapture 8.0×10⁻³ (~60× the old
flat-cap toy — 5 mm dome + ~21 mm neck/valve on-axis); **LS H-capture
1.2×10⁻² — the largest non-gas channel, 2.22 MeV γ born inside the
calorimeter** (not in any toy; scattered halo neutrons random-walk into the
four arms — see fig_scint_origin in the report). PlasticScint 2.2×10⁻³,
CFRP 3.5×10⁻⁴, Cu ~4×10⁻⁴ total.

---

## Data (all validated 2026-06-11: trees readable, exact entry counts)

`/eos/experiment/ntof/data/x17/full_sim/` (n_TOF experiment EOS — moved
from personal EOS after quota corruption, see below):

| dir | files | events | generator |
|---|---|---|---|
| `neutrons_subkev/` | 100 | 10⁹ | `--neutron`, E = 1 meV–1 keV |
| `neutrons_fullrange/` | 100 | 5×10⁸ | `--neutron`, E = 1 meV–100 MeV |
| `pairs_v2_step_target/` | 100 | 10⁷ | `--ipc 0.5` (X17+IPC 50/50), STEP geometry |

**Quota post-mortem:** the original 2026-06-10 campaign wrote to personal
EOS; the quota filled mid-campaign and **78 of 300 jobs produced GB-sized
ROOT files with no key directory** (look fine on disk, unreadable). All
were resubmitted with original seeds. **Lesson: after any campaign,
validate `uproot.open(f)['EventTree'].num_entries` on every file — `ls`
and byte counts lie.**

## Code & outputs (repo `MX17_Full_Geant`, all on main, lxplus clone at
`/afs/cern.ch/work/d/dneff/git/MX17_Full_Geant`)

- `scripts/analyze_thermal_captures.py` — EventTree scan: terminal-
  interaction budget, per-decade rates, direct (n,γ) list, absorption
  depths. Skips unreadable files. ~15 min for 10⁹ events with 6 workers.
  Output: `thermal_captures_subkev_full.npz/.json` (lxplus repo root,
  copies in `docs/report/`).
- `scripts/make_thermal_report_figures.py` — all report figures from the
  npz (+ `--scint-events <file>` for the two position-map figures).
- `scripts/plot_he3_xs.py` — ENDF cross-sections from `data/He3.h5`.
- `scripts/plot_geometry.py` — STEP-derived capsule (gas r=10 mm, 60 mm
  on-axis; Al 0.6 mm barrel / 5 mm dome / ~21 mm neck+valve; 0.9 mm CFRP),
  top-down + side view + 3D. Kept in sync with `DetectorConstruction.cc`.
- `docs/report/thermal_note.tex` → `thermal_note.pdf` (12 pp, final,
  full statistics). Build: pdflatex on lxplus (figures first; the He3.h5
  and EAR2 ROOT files live in `data/`).
- Condor submit scripts now default to the ntof EOS path.

## Next steps (where the next session should start)

1. **Talk to Alberto** (top priority, blocks ROI decision): which table
   rows anchored the sensitivity; IPC coefficient provenance; X17/IPC
   branching basis; capsule geometry update (table assumed 4 cm sphere).
2. **MeV-region capture budget** from run B-full (`neutrons_fullrange/`)
   — same scan machinery, wider energy binning; validates the thin-target
   rows where they claim to be valid.
3. If MeV ROI confirmed: **E_n-dependent pair kinematics** in
   `X17PrimaryGenerator`, re-plan run A statistics.
4. **Run C** (biased wall-background γs): `make_capture_library.py` on
   run-B output → `submit_neutrons.py --gamma-source capture_lib.csv`.
5. **Run D** (single-particle KE scans) + regenerate response JSON from
   run A (new STEP geometry) → rerun significance studies.
6. PLAN stages 3–4: `make_event_pools.py` slimming + Python pile-up
   sampler in `nTof_x17/MX17_Simulation`.
