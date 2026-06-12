# MeV-region campaign (active)

**Result (2026-06-12, run B-full, 5×10⁸ events, 1 meV–100 MeV):**
~33 X17 produced/day over the full range, **~22.5/day in the 0.2–2 MeV
window** (490 direct (n,γ) events, ±4.5% stat) — ×300 the sub-keV yield,
within 10% of the thin-target table row where the table is valid.
The measurement is feasible; remaining questions are acceptance and
backgrounds, not production.

**Task 2 (2026-06-12): pairs acceptance done.** X17 double-trigger 12.4%
(IPC 3.7%) at-rest baseline → **~2.8 double-triggered X17/day in-window**,
trigger-level S/B ≈ 0.09 before mass/angle cuts. Response JSON (10⁷ events)
installed in `nTof_x17/MX17_Simulation/`. Event pools (PLAN Stage 3) built
by `scripts/make_event_pools.py` → `analysis/pools/`.

| file | what |
|---|---|
| `HANDOFF_STAGE4_PILEUP.md` | **START HERE next session** — pile-up sampler task, pool format, in-flight builds |
| `HANDOFF_MEV_ANALYSIS.md` | the original handoff: task list, data, pitfalls, session log |
| `../../docs/report/mev_note.pdf` | the deliverable note (8 pp) |
| `../../analysis/mev/mev_rates.json` | per-decade rates, window integrals, errors |
| `../../analysis/mev/mev_captures.npz` | raw scan histograms + 714 (n,γ) energies |
| `../../analysis/pairs_v2/` | pairs PDF (17 sections) + response JSON |
| `../../analysis/pools/` | slimmed event pools for the Stage 4 pile-up sampler |

**Critical pitfall found this session:** the thermal
σ_nγ/σ_np = 1.0×10⁻⁸ must never be applied outside the sub-keV region —
the ratio rises to ~10⁻⁴ at MeV. Use the **direct (n,γ) counts** (the
fullrange run has plenty). An earlier pass that applied the thermal ratio
everywhere underreported the MeV yield by ~300×.

**Normalisation:** fullrange n/pulse = 2.263×10⁷ (integral of
`flux_n_pulse_NOisolet_100bpd`, 1 meV–100 MeV; sub-keV part reproduces the
7.31×10⁶ anchor). Alberto's table assumes 3.29×10⁷ — 31% hotter; compare
per incident neutron.
