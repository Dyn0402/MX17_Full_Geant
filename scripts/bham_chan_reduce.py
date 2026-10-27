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
q [MeVee]: the light, quenched per hit: e±/γ as deposited; protons by the Cecil et al. (NIM 161 (1979) 439)
plastic curve L = 0.95 E − 8.0 (1 − exp(−0.1 E^0.9)); deuterons/tritons the same at E/2, E/3 (scaled range);
α, ³He and heavier ions 0.02·E (a stand-in for their < 5 % light in plastic at < 2 MeV).
Legs, deposited and quenched: per arm, scintillator hits in time clusters (gap > 50 ns, as nfs_beam_reduce.py);
a leg is a SiPM bar ≥ 0.5 MIP (0.2376 MeV) AND a back plastic ≥ 0.5 MIP (1.675 MeV), in deposited energy (leg_e,
the trigger definition of every table) and in light (leg_q).  Rows lg_ev, lg_arm, lg_t, lg_e, lg_q.
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

acc = {k: [] for k in ("ev", "key", "t", "e", "q")}


def light(p, e):
    """per-hit electron-equivalent light [MeV] for particle names p and deposits e [MeV]."""
    q = e.copy()
    ion = ~np.isin(p, ["e-", "e+", "gamma", "mu-", "mu+"])
    for nm, k in (("proton", 1.0), ("deuteron", 2.0), ("triton", 3.0)):
        m = p == nm
        E = e[m] / k
        q[m] = k * np.clip(0.95 * E - 8.0 * (1 - np.exp(-0.1 * np.power(E, 0.9))), 0, None)
        ion &= ~m
    q[ion] = 0.02 * e[ion]
    return q
with uproot.open(a.file) as f:
    et = f["EventTree"]
    br = ["eventID", "neutron_E_eV"] + (["ev_weight"] if "ev_weight" in et.keys() else [])
    E = et.arrays(br, library="np")
    for t in f["HitTree"].iterate(["eventID", "armID", "detType", "particle", "u", "edep", "time"], library="np",
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
        pn = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in t["particle"][m]])
        acc["q"].append(light(pn, t["edep"][m] * 1e-6))
# legs per arm time cluster (deposited and quenched)
def legs(ev, key, tt, e, qq):
    arm, ch = key // 24, key % 24
    s = ch < 22                                          # SiPM bars + back plastics
    ev, arm, ch, tt, e, qq = ev[s], arm[s], ch[s], tt[s], e[s], qq[s]
    o = np.lexsort((tt, arm, ev))
    ev, arm, ch, tt, e, qq = ev[o], arm[o], ch[o], tt[o], e[o], qq[o]
    if not len(ev):
        return [np.zeros(0)] * 5
    nw = np.ones(len(ev), bool)
    nw[1:] = (ev[1:] != ev[:-1]) | (arm[1:] != arm[:-1]) | (np.diff(tt) > 50.0)
    c = np.cumsum(nw) - 1
    nc = c[-1] + 1
    out = []
    for x in (e, qq):
        bar = np.zeros((nc, 20))
        np.add.at(bar, (c[ch < 20], ch[ch < 20]), x[ch < 20])
        pl = np.bincount(c[ch >= 20], x[ch >= 20], minlength=nc)
        out.append((bar.max(1) >= 0.2376) & (pl >= 1.675))
    f = np.nonzero(nw)[0]
    k = out[0] | out[1]
    return ev[f][k], arm[f][k], tt[f][k], out[0][k], out[1][k]


ev, key, tt, e, qq = (np.concatenate(acc[k]) if acc[k] else np.zeros(0) for k in ("ev", "key", "t", "e", "q"))
lg = legs(ev, key, tt, e, qq)
o = np.lexsort((tt, key, ev))
ev, key, tt, e, qq = ev[o], key[o], tt[o], e[o], qq[o]
new = np.ones(len(ev), bool)
if len(ev):
    new[1:] = (ev[1:] != ev[:-1]) | (key[1:] != key[:-1]) | (np.diff(tt) > GAP_NS)
cid = np.cumsum(new) - 1
first = np.nonzero(new)[0]
E_cl = np.bincount(cid, e) if len(e) else e
Q_cl = np.bincount(cid, qq) if len(e) else e
ev_cl, key_cl, t_cl = ev[first], key[first], tt[first]
w_ev = np.ones(int(E["eventID"].max()) + 1) if len(E["eventID"]) else np.ones(1)
en_ev = np.zeros_like(w_ev)
if "ev_weight" in E:
    w_ev[E["eventID"]] = E["ev_weight"]
en_ev[E["eventID"]] = E["neutron_E_eV"] * 1e-6
np.savez_compressed(a.out, ev=ev_cl.astype(np.int32), arm=(key_cl // 24).astype(np.int8),
                    ch=(key_cl % 24).astype(np.int8), t=t_cl.astype(np.float64), e=E_cl.astype(np.float32), q=Q_cl.astype(np.float32),
                    w=w_ev[ev_cl].astype(np.float32), En=en_ev[ev_cl].astype(np.float32),
                    n_generated=np.array([len(E["eventID"])]),
                    lg_ev=lg[0].astype(np.int32), lg_arm=lg[1].astype(np.int8), lg_t=lg[2], lg_e=lg[3], lg_q=lg[4],
                    lg_w=w_ev[lg[0].astype(np.int64)].astype(np.float32) if len(lg[0]) else np.zeros(0, np.float32),
                    lg_En=en_ev[lg[0].astype(np.int64)].astype(np.float32) if len(lg[0]) else np.zeros(0, np.float32))
print(f"{a.file}: {len(E['eventID'])} events, {len(E_cl)} channel clusters -> {a.out}")
