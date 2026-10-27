#!/usr/bin/env python3
"""
bham_chan_reduce.py -- per-channel time clusters of one beam ROOT file, for the scintillator load
(x17_facility_search/birmingham/TARGET_GEANT.md; the beam analogue of lnl l10_singles.py).

    python3 bham_chan_reduce.py <neutrons_jobNNN_t0.root> -o <parts_chan/stem.npz>

Channels per arm (as l10_singles.py):
  ch 0-19  SiPM-wall bars (PlasticScint, 25 mm in u)
  ch 20    BackScintL (as built: the left 20×30×2 cm bar; big plastics: the 75×75×5 slab)
  ch 21    BackScintR (as built only)
  ch 22    LiqScint_1 (as built only)
  ch 23    Micromegas drift gas (DriftGas)
Each channel's hits are sorted in time; a gap > GAP_NS starts a new cluster (beam histories last µs to
s: moderated neutrons, capture γ, decays).  One row per cluster: event, arm, channel, t0 [ns], E [MeV],
and the event weight (EventTree ev_weight; 1 for analog runs), with the primary energy En [MeV].
"""
import argparse

import numpy as np
import uproot

GAP_NS = 100.0
CH = {"BackScintL": 20, "BackScintR": 21, "LiqScint_1": 22, "DriftGas": 23}

ap = argparse.ArgumentParser()
ap.add_argument("file")
ap.add_argument("-o", "--out", required=True)
a = ap.parse_args()

acc = {k: [] for k in ("ev", "key", "t", "e")}
with uproot.open(a.file) as f:
    et = f["EventTree"]
    br = ["eventID", "neutron_E_eV"] + (["ev_weight"] if "ev_weight" in et.keys() else [])
    E = et.arrays(br, library="np")
    for t in f["HitTree"].iterate(["eventID", "armID", "detType", "u", "edep", "time"], library="np",
                                  step_size="300 MB"):
        det = np.array([x.decode() if isinstance(x, bytes) else x for x in t["detType"]])
        ch = np.full(len(det), -1, np.int64)
        sip = det == "PlasticScint"
        ch[sip] = np.clip((t["u"][sip] + 250.0) // 25.0, 0, 19)
        for k, v in CH.items():
            ch[det == k] = v
        m = (ch >= 0) & (t["edep"] > 0)
        acc["ev"].append(t["eventID"][m].astype(np.int64))
        acc["key"].append(t["armID"][m].astype(np.int64) * 24 + ch[m])
        acc["t"].append(t["time"][m].astype(np.float64))
        acc["e"].append(t["edep"][m] * 1e-6)
ev, key, tt, e = (np.concatenate(acc[k]) if acc[k] else np.zeros(0) for k in ("ev", "key", "t", "e"))
o = np.lexsort((tt, key, ev))
ev, key, tt, e = ev[o], key[o], tt[o], e[o]
new = np.ones(len(ev), bool)
if len(ev):
    new[1:] = (ev[1:] != ev[:-1]) | (key[1:] != key[:-1]) | (np.diff(tt) > GAP_NS)
cid = np.cumsum(new) - 1
first = np.nonzero(new)[0]
E_cl = np.bincount(cid, e) if len(e) else e
ev_cl, key_cl, t_cl = ev[first], key[first], tt[first]
w_ev = np.ones(int(E["eventID"].max()) + 1) if len(E["eventID"]) else np.ones(1)
en_ev = np.zeros_like(w_ev)
if "ev_weight" in E:
    w_ev[E["eventID"]] = E["ev_weight"]
en_ev[E["eventID"]] = E["neutron_E_eV"] * 1e-6
np.savez_compressed(a.out, ev=ev_cl.astype(np.int32), arm=(key_cl // 24).astype(np.int8),
                    ch=(key_cl % 24).astype(np.int8), t=t_cl.astype(np.float64), e=E_cl.astype(np.float32),
                    w=w_ev[ev_cl].astype(np.float32), En=en_ev[ev_cl].astype(np.float32),
                    n_generated=np.array([len(E["eventID"])]))
print(f"{a.file}: {len(E['eventID'])} events, {len(E_cl)} channel clusters -> {a.out}")
