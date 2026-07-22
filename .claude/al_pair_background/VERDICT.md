# Al pair-production background — verdict

**2026-07-22 · full-statistics detector-level result.** Author: Claude (Opus 4.8) + Dylan.
Method + provenance in [PROGRESS.md]; figures in `analysis/al_pair/`.

## Bottom line

The ²⁷Al(n,γ) 7.72 MeV capture-γ pair-production background to the thermal X17
measurement is **enormous (6×10⁵ MM-reconstructed pairs/day) and its opening angle
is NOT a usable discriminant** — multiple scattering of the soft Al leptons drags
~7% of them past the X17 threshold, leaving Al outnumbering X17 by **~10⁷:1 inside
the X17 angular window**. Opening angle alone is hopeless.

**The background is defeatable by a tracking-based TOTAL pair-energy cut.** Al(7.72)
deposits ~6.7 MeV total lepton KE; X17 and IPC both carry the full 20.58 MeV
transition → ~19.6 MeV. A cut at ~13 MeV cleanly separates them. Requirement:
**per-lepton momentum resolution σ(p)/p ≲ 30%** from the MM (multiple-scattering
or range — there is no calorimeter). At that resolution the residual Al under the
X17 peak falls below the X17 signal rate while keeping ~87–95% of X17.

**So the thermal X17 measurement is feasible against Al IFF the MM delivers ≲30%
per-lepton momentum resolution.** That is the single make-or-break detector
requirement this study isolates. Note the softer-lepton-only cut is far weaker
(needs ≲10%); the total-energy cut is the one that works.

## Numbers (Deliverables 1–3)

Rates: in-gate (thermal) PRODUCTION × MM 2-track acceptance (both leptons make a
DriftGas hit). Al = 71,364 conv pairs from 10⁹ thermal neutrons; signal from
pairs_thermal_trig_2cm (300k events).

| | production/day | MM acc | MM pairs/day | reco θ>108° |
|---|---|---|---|---|
| **Al pair (capsule)** | 5.95×10⁶ | 10.9% | 6.5×10⁵ | 4.7×10⁴/day |
| IPC | 1.39 | 32.0% | 0.44 | 0.063/day |
| X17 | 0.035 | 33.0% | 0.012 | 0.0078/day |

- **Al truth opening angle** median 22°; **reco (MSC) median 46°**, tail >108° goes
  0.31% (truth) → **7.2% (reco)** — MSC is the whole story for the tail.
- **Soft-lepton kinematics:** Al softer lepton <3.35 MeV always (7.72 line shares
  6.70 MeV KE); X17 softer lepton >4 MeV always (massive boson). Clean gap.
- **S/B in the X17 window, opening angle only:** Al:X17 ≈ 6×10⁶ : 1.

## The discriminants (Deliverables 4–5)

| cut | needs | Al residual in X17 window | verdict |
|---|---|---|---|
| opening angle (θ>108°) | tracking dirs | 4.7×10⁴/day | hopeless (10⁷×X17) |
| + softer lepton >4 MeV | σ(p)≲10% | 49/day @10%, 950/day @20% | weak — thin 0.65 MeV gap |
| + **total pair E >13 MeV** | **σ(p)≲30%** | **<X17 for σ≲32%** (0.0017/day @30%) | **works** |

Analytic (Gaussian per-lepton resolution) total-energy leak: 1e-8/day @20%,
1.7e-3/day @30%, 0.2/day @40% → crosses X17 (0.0078/day) at **σ≈32%**.

## Caveats (must resolve before quoting a final significance)

1. **Momentum-resolution tails.** The 32% threshold assumes GAUSSIAN resolution.
   Real MSC/range momentum has non-Gaussian tails that push Al up — the true
   requirement needs the actual MM momentum-reconstruction performance, not a
   Gaussian model. This is the #1 follow-up.
2. **IPC is the irreducible background, not Al.** IPC total energy = 19.6 MeV too,
   so the total-energy cut does NOT remove it. IPC sits at **0.063/day in the X17
   window vs X17 0.0078/day (~8:1)** — that is the real signal-extraction problem,
   and its large-angle tail depends on the **IPC multipolarity, which the generator
   does NOT model** (1/Mee + isotropic decay, no E0/M1/E2). Reconcile with the E0
   branch before any significance. [[significance-projection]]
3. **Other capture lines.** N(n,γ) 10.8 MeV (total KE 9.78 MeV, 2.4×10⁴/day) is the
   closest line to the 13 MeV cut — still cleanly below it, but at the margin it
   needs slightly tighter resolution than the 7.72 line; 250× rarer than Al.
   All detector-material capture γ ≤ ~10.8 MeV → total KE ≤ ~9.8 MeV, all below 13.
4. **MC sample floor.** Direct MC proves Al suppression only to ~80/day (1 pair);
   the sub-X17 numbers are the analytic-tail extrapolation.
5. Same-arm fraction (Al 60%, X17 4%): forward Al pairs mostly enter one MM arm →
   2-track separation adds a further (unquantified here) Al inefficiency in our
   favour; wrong-track combinatorics (2000 SiPM/900 plastic singles/pulse) work
   against us — neither folded into the above.

## Feed-forward
- Get the real MM per-lepton momentum resolution (MSC fit / range) → replace the
  Gaussian model → firm σ requirement and residual Al.
- Fold residual Al + the IPC (multipolarity-correct) into the significance
  projection. Al is likely sub-dominant to IPC once the total-energy cut is in.
