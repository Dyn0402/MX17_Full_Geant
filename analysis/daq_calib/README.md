# DAQ calibration snapshots (pulled from the beam DAQ machine)

Exported from `daq_lxplus` (`~/daq/calibrations/` and
`~/PycharmProjects/nTof_x17/mx_july_beam_qa/calib/`) on 2026-07-19 so the
Geant4 trigger analysis is reproducible against the measured detector scales.
**These are snapshots — regenerate on the DAQ side after any HV/hardware
change; do not hand-edit.**

| file | what |
|---|---|
| `pss_mip_calib_run224489.json` | plastic per-PMT MIP peak (mV) + `mv_per_mev`, nominal HV, from WAL×PSS×LIQ triple-tagged spectra |
| `y88_energy_calib.json` | Y88 699/1612 keV Compton-edge energy scale (mV/keVee) per PMT/wall/LIQ channel (raised-HV source runs 224476–79) |
| `y88_vs_beam_mip.json` | Y88 edge vs beam-MIP cross-check (walls) |
| `mip_run224460.json` | SiPM-wall per-channel MIP peak (ADC & mV) |
| `threshold_scan_run224460.csv` | SiPM-wall trigger threshold → per-group MIP eff / purity / rate scan |
| `thresholds_run224460.json` | recommended SiPM-wall thresholds (12–14 mV, ≥95% weakest-group MIP eff) |

## Cross-check with the Geant4 sim (2026-07-19)

**MIP scale validated.** Sim plastic MIP = 4.33 MeV (Landau MPV of per-event
edep, 2.5 cm PVT); beam calibration uses mean dE = 5.046 MeV. Ratio
4.33/5.05 = 0.86 = the expected Landau MPV/mean for this thickness → the sim
reproduces the measured plastic energy scale. Sim SiPM-bar MIP (458 keV) is
consistent with 0.3 cm PVT. **Consequence: thresholds quoted in MIP transfer
directly between sim and hardware.**

**Threshold → hardware mV.** With the sim's optimum at plastic ≈ 0.6 MIP
(2.6 MeV), per-PMT thresholds = 0.6 × (measured MIP peak mV):

| PMT | MIP mV | 0.6 MIP mV | | PMT | MIP mV | 0.6 MIP mV |
|---|---|---|---|---|---|---|
| PSSA1 | 10.34 | 6.2 | | PSSC1 | 6.49 | 3.9 |
| PSSA2 | 8.06 | 4.8 | | PSSC2 | 6.32 | 3.8 |
| PSSB1 | 3.67 | 2.2 | | PSSD1 | 5.02 | 3.0 |
| PSSB2 | 5.94 | 3.6 | | PSSD2 | 3.27 | 2.0 |

**Discriminator floor ≈ 1.5 mV** (Y88 handoff) ⇒ the weakest PMT (PSSD2)
cannot set a threshold below ~0.46 MIP ≈ 2.0 MeV at nominal HV. This coincides
with the S/N optimum (~2.6 MeV), so there is no conflict; chasing lower
thresholds (to push singles to 300–400/window) requires raising plastic HV.

**Gain spread 3.2×** (PMT MIP peaks 3.3–10.3 mV) ⇒ set **per-PMT** mV, not one
global value, or a uniform mV gives a 3.2× spread in effective MIP threshold
across arms. Equalize first with `pss/global_gain_slide` if a common mV is
wanted.

**SiPM wall** runs 12–14 mV ≈ 0.25 of its MIP-sum peak (46–70 mV) at ≥95%
weakest-group MIP efficiency — well below the 0.5 MIP used in the sim scan, and
the sim shows signal is insensitive to the SiPM threshold over 0.1–0.5 MIP, so
no change to the recommendation.

## Open / to re-anchor

The 200/300/400-singles/window vertical lines in `../trigger_thermal/
threshold_linear.pdf` use the **simulated** singles rate (~450–550/window at
the floor). A clean scaler-based singles-vs-plastic-threshold measurement was
not yet available (the live `beam_july/test` plastic-threshold rate scan had
not written analysis products; M5 scaler dead in the July run). Re-anchor the
lines when that measurement lands.
