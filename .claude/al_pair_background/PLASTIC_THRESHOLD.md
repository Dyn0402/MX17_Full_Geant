# Plastic trigger-threshold optimization (thermal peak, DAQ-limited)

**2026-07-22.** Full Geant4 data: signal `pairs_thermal_trig_2cm` (2.5M ev),
background `neutrons_thermal_trig_2cm` (250M ev). Figure `analysis/al_pair/
plastic_opt.png`; scripts `digest_plastic_opt.py`, `plot_plastic_opt.py`.
MIP (2 cm, scaled ±20%): SiPM 0.458 MeV, **plastic 3.47 MeV**. Fixed SiPM leg at
0.5 MIP; scanned plastic. "Track" = charged DriftGas(MM) deposit ≥5 keV.

## Recommendation (updated 2026-07-22 for the 10 GbE DAQ upgrade)

**DAQ change:** the 10 GbE upgrade + readout start pushed to 1 ms raised the
thermal-band readout budget from ~4 to **~24 events/spill** (`nTof_x17_DAQ`
`docs/network_upgrade_10g/04_bandwidth_model.md`, IPD 10; ~20 at the IPD-15
fallback). The 1 ms start matches the sim's t>1 ms gate, so background rates are
unchanged — only the budget moved.

**Raise the plastic threshold from 0.5 MIP to ~1.4 MIP (~4.9 MeV)** — set it so
the in-gate background trigger rate equals the new ~24/pulse budget. Expected:
**≈4.8× more IPC recorded** than running 0.5 MIP at the same budget (and ~28× vs
the old 0.5-MIP / 4-per-pulse operating point). IPC efficiency at the optimum is
4.5% (vs 2.5% at the old, tighter 1.75-MIP recommendation) — the bigger budget
buys back efficiency by letting you sit at a lower threshold, right at the IPC
plastic-MIP turn-on.

Robust operational rule (independent of the ±20% MIP calibration and the exact
budget): **tune the plastic threshold on the measured thermal-gate trigger rate
until it equals the DAQ budget.** The optimum sits exactly there.

**Optimum vs budget** (bg = 202/pulse at 0.5 MIP; IPC-read gain vs 0.5 MIP at
that budget):

| DAQ budget/pulse | opt plastic | IPC eff (≥1 leg) | IPC recorded gain |
|---|---|---|---|
| 4 (old) | 1.76 MIP (6.1 MeV) | 2.5% | ×16 |
| 12 | 1.57 MIP (5.5 MeV) | 3.4% | ×7.4 |
| **24 (current best)** | **1.41 MIP (4.9 MeV)** | **4.5%** | **×4.8** |
| 48 | 1.23 MIP (4.3 MeV) | 6.1% | ×3.3 |
| 96 | 0.98 MIP (3.4 MeV) | 7.6% | ×2.0 |

As the budget grows the optimum falls toward ~1 MIP and the threshold-tuning gain
shrinks (you're less saturated). If the switch only negotiated 2.5 GbE (budget
~11–15), read off the 12-budget row (~1.55 MIP).

## Why (the counter-intuitive part)

IPC is **~10⁻⁴ triggers/pulse** — vanishingly rare against the **~200/pulse**
background trigger stream. At the current 0.5 MIP the DAQ (4/pulse) is saturated
~50×, so it reads a random 4 and **captures only ~2% of the IPC that triggers**.
So the goal is NOT to enrich IPC in the 4 reads (it is far too rare) — it is to
**cut the background trigger rate down to the DAQ budget so the DAQ is no longer
saturated when a real IPC fires.**

| plastic thr | IPC eff (≥1 leg) | bg trig/pulse | DAQ read-frac | IPC recorded |
|---|---|---|---|---|
| 0.5 MIP (1.73 MeV, current) | 7.8% | 202 | 2% | ×1 (baseline) |
| **1.76 MIP (6.1 MeV, optimum)** | 2.5% | **4.0** | **100%** | **×16** |

The optimum is exactly where bg = budget: below it the DAQ saturates and drops
IPC randomly; above it you sacrifice IPC efficiency for no throughput gain.
Raising the threshold costs raw IPC efficiency (7.8→2.5%) but the ×50 gain in
read-fraction wins by ×16 net.

The physics that makes a high plastic cut keep IPC: **IPC/X17 leptons are
energetic (up to 19 MeV) and punch through the plastic as a MIP (~4 MeV peak,
Panel A), while the Al backgrounds — both double-Compton electrons AND the soft
Al-pair leptons (<3.35 MeV) — deposit little in the plastic.** So a plastic
threshold near/above the MIP preferentially keeps the energetic signal.

## Important caveat on "minimize double-Compton"

The plastic threshold does **not selectively remove double-Compton.** Both
double-Compton and real tracks put an *electron* in the plastic, so their plastic
spectra overlap; double-Compton is still ~85% of the 4 reads even at the optimum.
What the threshold does is cut the **total** rate to de-saturate the DAQ.
**Selective double-Compton rejection requires the MM track — which is not
available at the hardware trigger** (the TPC is too slow). So double-Compton is
removed OFFLINE (require an MM track), not by the plastic threshold; the trigger's
only job is to not saturate. (Downstream, the total-energy cut from [VERDICT.md]
then separates IPC from the surviving Al-pair tracks.)

## Caveats
- MIP calibration is the 2.5→2.0 cm scaled value (±20%). Quote the threshold in
  MIP / as a rate-based rule, not a fixed MeV. A real 2 cm muon calibration would
  firm up the MeV value.
- DAQ model = random sampling of triggers under saturation (read-frac =
  min(1, budget/rate)); reasonable since IPC and Al share the thermal TOF time
  profile, but the true first-N-then-deadtime behaviour should be checked.
- Assumes the trigger = ≥1 leg (any arm). A 2-arm coincidence trigger lowers the
  background rate (but double-Compton can still fake 2 arms via cross-arm scatter)
  and would shift the optimum lower — worth a follow-up scan.
- The DAQ budget (4/pulse) drives the optimum: a larger budget → lower threshold →
  higher IPC efficiency. Panel C peak moves with the budget.
