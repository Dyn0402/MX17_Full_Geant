#!/usr/bin/env python3
"""
bham_chan_merge.py -- bham_chan_reduce.py parts of one sample -> weighted per-channel spectra
(x17_facility_search/birmingham/TARGET_GEANT.md §5, the scintillator load).

    python3 bham_chan_merge.py <parts_chan/*.npz> -o <sample>_chan.npz

Rows are channel time clusters (gap > 100 ns) with the event weight w and the primary energy En.
Output, all Σ w (divide by n_generated for per-primary rates):
  h[iEn, arm, ch, iE]    clusters per primary-energy bin, arm, channel, deposited-energy bin
  p[iEn, arm, ch]        Σ w·E deposited [MeV]
  hq, pq                 the same in quenched light (MeVee; parts from the quenching reducer only)
  en_edges (MeV, 0-20 in 0.05 below 1.3, then coarse), e_edges (log, 1 keV-200 MeV); n_generated.
"""
import argparse

import numpy as np

EN = np.unique(np.r_[np.arange(0.0, 1.3001, 0.05), [2, 5, 10, 14, 15, 16, 17, 18, 19, 20.5]])
EE = np.geomspace(1e-3, 200.0, 85)

ap = argparse.ArgumentParser()
ap.add_argument("parts", nargs="+")
ap.add_argument("-o", "--out", required=True)
a = ap.parse_args()
h = np.zeros((len(EN) - 1, 4, 24, len(EE) - 1))
p = np.zeros((len(EN) - 1, 4, 24))
hq, pq = np.zeros_like(h), np.zeros_like(p)
ngen = 0
for f in a.parts:
    z = np.load(f)
    ngen += int(z["n_generated"][0])
    if not len(z["e"]):
        continue
    ie = np.clip(np.digitize(z["En"], EN) - 1, 0, len(EN) - 2)
    ik = np.clip(np.digitize(z["e"], EE) - 1, 0, len(EE) - 2)
    np.add.at(h, (ie, z["arm"].astype(int), z["ch"].astype(int), ik), z["w"])
    np.add.at(p, (ie, z["arm"].astype(int), z["ch"].astype(int)), z["w"] * z["e"])
    if "q" in z.files:
        iq = np.clip(np.digitize(z["q"], EE) - 1, 0, len(EE) - 2)
        np.add.at(hq, (ie, z["arm"].astype(int), z["ch"].astype(int), iq), z["w"] * (z["q"] > EE[0]))
        np.add.at(pq, (ie, z["arm"].astype(int), z["ch"].astype(int)), z["w"] * z["q"])
np.savez_compressed(a.out, h=h, p=p, hq=hq, pq=pq, en_edges=EN, e_edges=EE, n_generated=np.array([ngen]))
print(f"{len(a.parts)} parts, {ngen} primaries -> {a.out}")
