# Handoff: thermal-gate (>1 ms) trigger optimization — 2026-07-19

Session goal (Dylan, 2026-07-18): with the final surveyed geometry, determine
what the detectors see in the thermal-neutron time range (>1 ms after the γ
flash) and optimize the trigger between the **SiPM wall** and the **plastics**
to tag IPC/X17-like pairs (self-shielding included) while rejecting capture-γ
backgrounds. Status: **complete** — results in `CAMPAIGN_STATUS.md`
(2026-07-19 headline) and below. One bookkeeping pass still running (§6).

## 1. The answer

**Recommended menu: 2 SiPM-wall legs ≥ 0.5 MIP (per-bar discriminator) +
≥ 1 plastic bar ≥ 1.0–1.2 MIP.**
ε(X17) = 7.0→4.0 %, ε(IPC | both legs > 5 MeV) = 5.3→3.2 %, background
1.8→0.8 pair-tags/pulse in the gate (≈0.4→0.2 Hz — DAQ-trivial).
- MIP scales (sim MPV, normal-incidence µ⁻): **SiPM bar 458 keV, plastic
  4.33 MeV**. Dylan will calibrate hardware the same way (MIP peaks), so
  thresholds transfer directly in MIP units.
- Full SiPM×plastic coincidence on both legs costs ~5× signal (plastics
  solid angle) → use only if purity trumps everything: (0.5, 0.5–0.8) MIP
  gives ε 2.1–1.3 %, bg 0.17–0.03/pulse.
- Plastic threshold physics: signal legs punch through 2.5 cm PVT at ~1 MIP
  (stopping legs up to ~3 MIP); H-capture 2.22 MeV γs < 0.5 MIP; thermal bg
  ends at ~1.7 MIP (7.72 MeV ²⁸Al Compton edge). **Do not go ≥ 1.5 MIP** —
  ε collapses to 1.4 %.
- Epithermal spill-in (delayed H-captures drifting into the gate):
  negligible — 0 pair-tags / 5×10⁸ epi neutrons, singles ~8/pulse vs ~170
  thermal.
- Caveat: in-gate ³He(n,γ) sits on the self-shielding ceiling
  (~1×10⁻⁴ IPC/pulse) — the >1 ms gate is background/veto territory; the
  physics ROI remains the MeV window (see 2026-06-12 headline).

## 2. Data (all validated: tree readable + exact entry counts)

On `/eos/experiment/ntof/data/x17/full_sim/`:

| Dir | Contents | Norm |
|---|---|---|
| `neutrons_thermal_trig/` | 100×10⁷ EAR2-flux neutrons, 1 meV–2 eV | 4.284×10⁶ n/pulse in window |
| `neutrons_epi_trig/` | 50×10⁷ neutrons, 2 eV–100 keV | 6.484×10⁶ n/pulse |
| `pairs_thermal_trig/` | 100×10⁵ X17+IPC (50/50), self-shielded vertices | — |

Gate ⇔ energy: TOF[ms] = 1.41/√E_n[eV] at 19.5 m; >1 ms ⇔ E_n < 1.99 eV.
Flux integrals from `flux_n_pulse_NOisolet_100bpd` in
`data/fluxEAR2-Ph3_in_different_units.root` (full range 2.263×10⁷ n/pulse).

On lxplus (`/afs/cern.ch/work/d/dneff/git/MX17_Full_Geant`):
- `analysis/trigger_thermal/`: `trigger_scan.json` (the full 20×20 scan),
  `trigger_spectra.pdf`, `trigger_roc.pdf`, `run.log` (headline table),
  `capture_lib_thermal.csv` + `gas_vtx_lib_full.csv` (pending, §6).
- `smoke/`: `mip_mu_armD_t0.root` (MIP calibration, 20k µ⁻ 1 GeV, θ=90°
  φ=3° — φ≠0 avoids the bar-9/10 boundary), `thermal_gas_vtx_lib.csv`
  (the 155k-vertex library used by the pairs campaign, from a 200k-neutron
  local run), assorted smoke outputs.

Locally: `analysis/trigger_thermal/` holds copies of the scan JSON + PDFs.

## 3. New machinery (commits 808ffe0..HEAD)

- **`--pair-vertex-lib <csv>`** (`mx17_full_sim`): X17/IPC vertices sampled
  from He3Gas capture positions instead of uniform-in-gas. This is how
  thermal self-shielding enters the signal: measured captures sit at median
  14 mm depth, 95 % within 25 mm of the gas entrance — the pair source is
  ~3.7 cm upstream of target centre (visible acceptance asymmetry).
- **`make_capture_library.py --gas-lib <csv> --gas-emax <eV>`**: writes the
  He3Gas capture-position library ((n,p) positions trace the same attenuated
  spatial law as the rare (n,γ)). Default emax 2 eV = the gate.
- **`submit_pairs.py --vertex-lib <csv>`**: pass-through.
- **`scripts/analyze_trigger_thermal.py`**: streams HitTrees (chunked, ok
  for 400 GB), builds per-(event, arm) observables S = max SiPM-bar edep,
  P = max plastic-bar edep (per-channel, MIP units), applies the per-hit
  TOF gate for epi files, scans A×B ∈ [0.1, 2.0]² MIP over five topologies
  (leg1, leg2, leg2_confirm1, SiPM-only, plastic-only). SiPM bar index from
  arm-local u: `floor((u_mm+250)/25)`. MIP MPV from `--mip` file(s).

Detector name → class map: `PlasticScint*` = SiPM wall bars (3 mm),
`BackScintL/R` = plastics (25 mm), `LiqScint_1` = LS (recorded, not in the
trigger logic — an open item, §7).

## 4. Reproduce / extend

```bash
# local → lxplus deploy (lxplus repo tracks the same GitHub origin now)
git push origin main            # local
ssh lxplus 'cd /afs/cern.ch/work/d/dneff/git/MX17_Full_Geant && git pull'
# build (lxplus only — no local Geant4)
source scripts/setup_lxplus.sh && cd build && cmake .. && make -j4

# the three campaigns (as submitted 2026-07-18)
python3 scripts/submit_neutrons.py --emin 1e-3 --emax 2   --njobs 100 --nevents 10000000 \
   --outdir .../neutrons_thermal_trig --jobdir ~/condor/mx17_neutrons_thermal_trig
python3 scripts/submit_neutrons.py --emin 2    --emax 1e5 --njobs 50  --nevents 10000000 \
   --outdir .../neutrons_epi_trig     --jobdir ~/condor/mx17_neutrons_epi_trig
python3 scripts/submit_pairs.py --njobs 100 --nevents 100000 --ipc 0.5 \
   --vertex-lib smoke/thermal_gas_vtx_lib.csv --outdir .../pairs_thermal_trig \
   --jobdir ~/condor/mx17_pairs_thermal_trig

# MIP calibration + analysis
./build/mx17_full_sim -n 20000 -t 1 -s 21 -o smoke/mip_mu_armD --single mu- 1000 90 3
python3 scripts/analyze_trigger_thermal.py \
   --signal .../pairs_thermal_trig --thermal .../neutrons_thermal_trig \
   --epi .../neutrons_epi_trig --mip smoke/mip_mu_armD_t0.root \
   -o analysis/trigger_thermal
```

## 5. Environment gotchas (will bite you)

1. **Run the sim single-threaded** (`-t 1`, and Condor jobs default to it).
   `-t ≥ 2` races ROOT's plugin manager in `RunAction::BeginOfRunAction`
   (concurrent per-thread `TFile::Open`) → worker crash. Fix would be
   `ROOT::EnableThreadSafety()`; not needed for production.
2. Output files are named `<base>_t0.root` (MT build, 1 worker). Globs
   handle it; exact-name code won't.
3. **Kerberos**: TGT lasts 24 h. Expiry silently wedges nohup'd lxplus jobs
   (AFS *and* EOS-fuse lost; processes survive but do nothing). Check
   `klist` before multi-hour runs. After `kinit dneff@CERN.CH`, restart the
   ssh master (`ssh -O exit lxplus`) — and note lxplus now asks a **2FA OTP
   on new ssh masters**, which needs Dylan at the keyboard. Condor jobs are
   immune. Wedged processes must be killed and restarted.
4. Validate every campaign file (EventTree.num_entries == expected) before
   analysis — `ls` is not enough (June quota incident).
5. lxplus `/tmp` is node-local; use the AFS work dir for anything another
   session must see.

## 6. Still running at handoff time

`make_capture_library.py` over the full 10⁹ thermal campaign (started
2026-07-19 02:11, EventTree-only; slow because of 10⁹ string branches):
writes `analysis/trigger_thermal/capture_lib_thermal.csv` (wall library),
`gas_vtx_lib_full.csv` (definitive 500k-row vertex library) and prints the
final per-volume capture budget to `capture_budget.log` **only at the end**
(empty log = normal). Check:
`ssh lxplus 'pgrep -f make_capture_library || cat .../capture_budget.log'`.
Expected (200k-n preview): He3 (n,p) 76.7 %, Al dome 6.5×10⁻³/n, LS H-capture
2.1×10⁻³/n, plastics ~1.4×10⁻³/n, SiPM bars 4×10⁻⁴/n. Nothing downstream
blocks on it.

## 7. Natural next steps

1. **MeV-window trigger scan** — same machinery, µs-scale gate, rates ×10⁴:
   this is the physics ROI. Needs a gate definition (1.0–3.2 µs at 19.5 m)
   and E_n-dependent pair kinematics (E* ≈ 20.58 + ¾E_n MeV) for the signal
   run — the generator extension is still the open item from
   PLAN_NEUTRON_CAMPAIGN.md step 7.
2. Add the LS to the trigger logic (currently recorded but unused) — e.g.
   as an additional confirm layer or offline energy sum.
3. Accidentals/pile-up: feed `trigger_scan.json` leg-singles rates (~170
   legs/pulse in-gate at loose b) into the Python fast-MC coincidence
   machinery (`MX17_Simulation`), which owns the time domain.
4. Refresh the response JSON / event pools for the fast-MC from
   `pairs_thermal_trig` (final geometry) — `make_event_pools.py` was the
   old Stage-3 plan.
5. When the full capture library lands, re-check the vertex-library shape
   against the 155k bootstrap (expected identical within stats).
