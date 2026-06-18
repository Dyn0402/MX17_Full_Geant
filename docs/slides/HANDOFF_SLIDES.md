# Handoff: results slide deck + significance projection (start here)

**Written:** 2026-06-18.
**For:** a future session returning to `docs/slides/` — the INFN/collaboration
deck, the rebinned spectra, and the (caveated) confidence-level projection.

---

## What this session did

A colleague asked for two things off the slide-4/5 opening-angle spectra:
(1) a **factor-4 rebin** of those spectra (kept in backup slides), and
(2) a **confidence level** for the X17 signal. Both are done; the CL number is
deliberately framed as a *statistical ceiling*, not a final result.

### 1. Rebinned spectra (×4, 2° → 8° bins)
- `scripts/make_slides_figures.py` gained a `rebin` parameter
  (`stacked_panel`, `fig_stacked`, `fig_compare`; helper `rebinned_edges`,
  const `BASE_BIN_DEG = 2.0`). `main()` now also emits the ×4 versions.
- New figures (originals **untouched**): `docs/slides/figs/`
  `fig_stacked_july_rebin4.png`, `fig_stacked_ls3_rebin4.png`,
  `fig_stacked_compare_rebin4.png`.
- Edge effect (documented in code): 90 base bins isn't divisible by 4, so the
  last bin (176–180°) holds the 2 leftover base bins. Distribution is empty
  there — cosmetic only.

### 2. Significance projection — deck slide 12
Computed with an **Asimov binned profile-likelihood ratio**:
`Z = sqrt(2·Σ[(s+b)ln(1+s/b) − s])` over the smeared opening-angle templates.
Quick calc lived in `/tmp/cl_calc.py` (not committed — see "redo it" below).

Results (30-day run; α_IPC=3.5e-3, X17/IPC=2.5%, MM-double acceptance,
best-estimator capsule-MS smearing):

| Scenario | S / B | Shape-fit Z (smeared) | Z (no smear) | naïve S/√B | best θ-window cut |
|---|---|---|---|---|---|
| July 0.2–0.7 MeV | 64 / 3093 | **2.6σ** | 3.9σ | 1.2 | 2.5σ ([98,158°]) |
| Post-LS3 0.2–2 MeV | 220 / 10600 | **4.9σ** | 7.3σ | 2.1 | 4.7σ ([98,158°]) |

Key findings:
- Smearing costs ~0.67× in Z; does **not** kill the measurement.
- **Binning-independent**: 2° vs 8° agree to <1% (so the rebin is purely
  cosmetic — good reassurance for the colleague).
- The broad shoulder makes a simple cut-and-count nearly as good as the full
  fit (within ~5%), *if* the background is known.
- Shape fit ~doubles the naïve whole-spectrum S/√B by weighting the high-θ
  shoulder where IPC is sparse.

### Deck assembly
- `scripts/build_slides.py` now builds **16 slides** (was 12). Added:
  - helpers `_redbox` (a KEEP-this red caveat box, distinct from the gold
    delete-before-INFN `_caveat`) and `_table`.
  - **Slide 12**: significance table + "What IS in these numbers" box +
    prominent red caveat box.
  - **Slides 14–16**: "Backup: factor-4 rebinned spectra" divider + rebinned
    July + rebinned July-vs-LS3.
- Rebuild: `.venv/bin/python scripts/build_slides.py` → overwrites
  `docs/slides/MX17_results.pptx`.

---

## THE CAVEAT (most important thing to remember)

**2.6σ / 4.9σ are a statistics-only ceiling. The IPC angular SHAPE is assumed
equal to the simulation truth. Its uncertainty will dominate the real CL and is
NOT propagated.** User (D.N.) was explicit that quoting a CL without these is
misleading — the deck slide 12 red box states this, but do not let the number
escape the caveat elsewhere. Not included:
1. **Physics shape** — single E1 transition assumed (probably fair, unquantified).
2. **Detector / acceptance shape distortion** — poor efficiency separating
   near-collinear tracks (small θ), and finite active-area acceptance falling
   off for wide-angle pairs (large θ). Both reshape the **high-θ region where
   the signal lives**.

Quantitative kicker: with B ≈ 3k–10k, a few-% shape error in the signal region
already rivals the statistical fluctuation → **by post-LS3 the IPC shape, not
statistics, sets the achievable CL.** Background *normalization* is cheap
(float it from the X17-free θ≲90° region, ~2%, negligible) — it's the shape,
not the norm, that matters.

---

## Where we left off / next step

The proper analysis (the real deliverable behind slide 12, and the
"template-fit significance projection" open item on the summary slide) is **not
done**:

- Build a toy-MC / profile-likelihood fit with an **IPC-shape nuisance
  parameter** + **full detector response** (track-separation efficiency vs θ,
  active-area geometric acceptance vs θ) folded into the IPC template.
- That turns the ceiling into a defensible median-expected Z and p-value vs run
  time. The `/tmp/cl_calc.py` Asimov calc is the skeleton to promote into a
  committed script (suggest `scripts/significance_projection.py`).

To regenerate everything from scratch:
```
.venv = /home/dylan/PycharmProjects/nTof_x17/.venv/bin/python
$VENV scripts/make_slides_figures.py    # figures incl. rebin4
$VENV scripts/build_slides.py           # -> docs/slides/MX17_results.pptx (16 slides)
```

**Note:** LibreOffice is broken in this environment (`libreglo.so` missing), so
the deck could not be rendered to image to eyeball slide 12's layout. It's
sized within bounds, but verify spacing in PowerPoint/Google Slides.

Inputs: `analysis/mev/mev_captures.npz`, `analysis/mev/mev_rates.json`,
`analysis/pairs_v2/geant4_response.json`. Templates come from
`scripts/make_angular_resolution_figs.py` (sample_pairs / smear_pair_response).

Memory: `significance-projection.md` in the auto-memory dir mirrors the headline
numbers + caveat.
