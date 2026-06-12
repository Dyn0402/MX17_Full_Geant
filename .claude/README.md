# .claude/ — session docs & handoffs

Technical markdown for the analysis campaigns, organized by phase. Plain-text
campaign-level docs (`CAMPAIGN_STATUS.md`, `PLAN_NEUTRON_CAMPAIGN.md`,
`HANDOFF_FULL_SIM.md`) stay at the repo root; pitfall notes referenced by the
reports stay in `docs/`.

| where | what |
|---|---|
| `thermal/` | sub-keV campaign handoff (concluded: measurement infeasible) |
| `mev/` | MeV-region campaign handoff + session notes (active) |
| `../analysis/<phase>/` | scan outputs (npz raw counts, derived rate JSONs) |
| `../docs/report/` | the deliverable notes: `thermal_note.pdf`, `mev_note.pdf` |
| `../scripts/` | `analyze_*_captures.py` (scans), `make_*_report_figures.py` (figs+rates) |
