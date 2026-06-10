# EAR2 Neutron Beam Data — n_TOF Phase 3

This directory contains the evaluated neutron flux and spatial beam profile data for
EAR2 (Experimental Area 2) at the n_TOF facility, CERN. These files are used to
drive the neutron beam input for the MX17 Geant4 simulation.

All files were retrieved from CERN EOS/AFS via lxplus (June 2025).

---

## ROOT Files

### `fluxEAR2-Ph3_in_different_units.root`
**The primary flux file.** Contains the Phase 3 evaluated EAR2 neutron flux in four
unit variants, all as `TH1D` with 1100 log-spaced bins from **1 meV to 100 MeV**
(100 bins/decade).

Source: `/eos/experiment/ntof/flux_RF_BIF/FLUX/EAR2/`  
Reference: Sabaté-Gilarte et al., evaluated flux for n_TOF Phase 3.

| Histogram name | Units | Notes |
|---|---|---|
| `flux_n_pulse_isolet_100bpd` | n / eV / pulse | Isolethargic (E·dΦ/dE); use for reaction rate integrals |
| `flux_n_pulse_NOisolet_100bpd` | n / pulse | Flux per bin; use for counting neutrons in an energy interval |
| `flux_n_cm2_pulse_100bpd` | n / cm² / pulse | Assumes uniform beam within r = 2.5 cm collimator |
| `flux_n_cm2_s_100bpd` | n / cm² / s | Same, converted at 0.8 Hz rep rate |

**Normalisation conventions:**
- All histograms assume **7×10¹² protons per pulse** on the spallation target.
- The cm²-normalised histograms assume a **uniform beam within r = 2.5 cm** — this is
  a simplification; use `lamda2DvsEn_EAR2.root` for the true spatial profile.
- ~3.2×10⁵ n/pulse pass through the full 2.5 cm collimator in the X17 ROI (10–20 MeV).

**For Geant4:** Read `flux_n_pulse_NOisolet_100bpd` to get the number of neutrons per
pulse per energy bin. Use this as a sampling distribution (or weighting function) for
the primary neutron energy. Combined with the spatial profile below, this fully
characterises the beam.

---

### `lamda2DvsEn_EAR2.root`
**The spatial beam profile file.** Contains one object:

- `Lambda2D` — `TH2D`, shape **(240 energy bins) × (3000 radial bins)**
  - X axis: neutron energy (eV), log-spaced, **1 meV to ~1 GeV**, 240 bins
  - Y axis: radial distance from beam axis (cm), **0 to 3.0 cm**, bin width 10 µm

Source: `/afs/cern.ch/work/m/msabateg/public/n_TOF-files/` (Marta Sabaté-Gilarte)  
Produced by: the n_TOF optical transport code (`ntof/transport`) run on FLUKA
spallation-target output, propagated to L = 19.5 m (EAR2 sample position).

**How to interpret:** Each row (energy bin) gives the radial distribution of neutrons
at the EAR2 sample position for that energy. The values are Monte Carlo weights
(not absolute flux); normalise each energy row by its sum to get the probability
density dP/dr. The radial bins are fine (10 µm), so rebin by ~50× before use
(0.5 mm resolution is adequate).

**Key results from this file:**
- At thermal/eV energies: beam is narrow, ~90% within r < 0.5 cm.
- At 1 MeV: beam has broadened, ~90% within r < 1.5 cm.
- In X17 ROI (10–20 MeV): ~47% within r < 0.5 cm, ~91% within r < 1.0 cm,
  essentially 100% within r < 1.5 cm.
- Above ~100 MeV: beam spreads to fill the full 3 cm scoring window.

**For Geant4:** To sample the (x, y) position of a primary neutron at a given energy E:
1. Find the energy bin index for E in the X axis.
2. Extract that row as the radial probability distribution dP/dr.
3. Sample r from this distribution (using inverse CDF or rejection sampling).
4. Sample φ uniformly in [0, 2π); convert to (x, y) = (r·cos φ, r·sin φ).

---

### `nTOF-Ph3_fluence_EAR2_2014-2015.root`
Older format of the same Ph3 evaluated flux. Contains one histogram:
`hFlux_eval_ear2_2014_2015_100bpd` — identical data to `flux_n_pulse_isolet_100bpd`
in the main file, stored as `TH1D` (TH1F format). Provided for cross-checks only;
use `fluxEAR2-Ph3_in_different_units.root` instead.

---

### `nTOF-Ph3_fluence_EAR2_2014-2015_TOF-E.root`
High-resolution (1000 bins/decade) version of the flux and its TOF-space equivalent.
Useful if you need finer energy binning than 100 bpd, e.g. near resonance structure.

| Histogram | Description |
|---|---|
| `hE_eval_ear2_2014_2015_1000bpd` | Isolethargic flux vs true energy, 11000 bins, 1 meV–100 MeV |
| `hTOF_SiMon2pos_ear2_2014_2015_1000bpd` | Flux binned in TOF space (SiMon2 detector), 5500 bins |

The TOF-space histogram has a distorted shape at high energies due to the Jacobian
of the TOF→E transformation — this is expected and not a calibration error.

---

### `fluxEAR2_TH1F.root`
Identical to `nTOF-Ph3_fluence_EAR2_2014-2015.root`, stored as `TH1F` rather than
`TH1D`. Kept for completeness.

---

### `fluka_flux_ear2.root`
FLUKA-simulated 1D EAR2 flux (`simulated_flux_EAR2`, `TH1D`). Useful for comparing
the evaluated flux against the pure simulation before experimental normalisation.

---

## Plotting Scripts

Both scripts require Python 3 with `numpy`, `matplotlib`, `scipy`, and `uproot`.
Run from within the `data/` directory using the nTof_x17 virtual environment:

```bash
python plot_ear2_flux.py
python plot_ear2_spatial_profile.py
```

Or with the explicit venv path:

```bash
~/PycharmProjects/nTof_x17/venv/bin/python plot_ear2_flux.py
```

### `plot_ear2_flux.py`
Produces four figures (saved as PNGs in this directory):

| Output file | Contents |
|---|---|
| `ear2_flux_units.png` | All four flux unit variants on log-log axes |
| `ear2_flux_ph3_vs_old.png` | Ph3 evaluated vs 2014–2015 with ratio panel |
| `ear2_flux_tof.png` | High-resolution (1000 bpd) and TOF-space flux |
| `ear2_flux_radial_context.png` | Cumulative flux above threshold; flux vs target radius (uniform beam) |

### `plot_ear2_spatial_profile.py`
Produces four figures from `lamda2DvsEn_EAR2.root`:

| Output file | Contents |
|---|---|
| `ear2_spatial_2d_map.png` | 2D colour map of radial profile vs energy (normalised per slice) |
| `ear2_spatial_radial_slices.png` | Radial profiles at 7 representative energies |
| `ear2_spatial_beam_size.png` | Mean r, RMS r, r₅₀, r₉₀ vs energy |
| `ear2_spatial_target_fraction.png` | Fraction of beam within r < r_target vs energy |

---

## Origin of Data on CERN Systems

| File | Source path on lxplus |
|---|---|
| `fluxEAR2-Ph3_in_different_units.root` | `/eos/experiment/ntof/flux_RF_BIF/FLUX/EAR2/` |
| `fluxEAR2_TH1F.root` | `/eos/experiment/ntof/flux_RF_BIF/FLUX/EAR2/Additional_files/` |
| `nTOF-Ph3_fluence_EAR2_2014-2015.root` | `/eos/experiment/ntof/flux_RF_BIF/FLUX/EAR2/Additional_files/` |
| `nTOF-Ph3_fluence_EAR2_2014-2015_TOF-E.root` | `/eos/experiment/ntof/flux_RF_BIF/FLUX/EAR2/Additional_files/` |
| `lamda2DvsEn_EAR2.root` | `/afs/cern.ch/work/m/msabateg/public/n_TOF-files/` |
| `fluka_flux_ear2.root` | `/afs/cern.ch/work/m/msabateg/public/n_TOF-files/` |

The spatial profile was produced by the n_TOF transport code
(`/eos/experiment/ntof/simul/transport/`) using FLUKA spallation-target phase-space
files from Phase 3, propagated to L = 19.5 m with the EAR2 collimator geometry.
The transport code executable and input files are at
`/afs/cern.ch/work/m/msabateg/public/transport_tools/transport_lxplus/`.
