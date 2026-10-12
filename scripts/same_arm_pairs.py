#!/usr/bin/env python3
"""
same_arm_pairs.py -- every lepton that enters a Micromegas drift gap, with its
ancestry, for the same-chamber pair study of nTof_x17
(sept26_prelim_analysis/same_chamber_physics.py, SAME_CHAMBER_PAIRS.md).

The question it serves: of the events with TWO leptons in one chamber's gap,
how many are an e+e- pair (one converting gamma, or an internal pair), how many
are two Compton electrons of two gammas of one capture, how many a Compton
electron and its own delta ray -- and how many come with a trigger.

    # per ROOT file (lxplus / condor, LCG_106)
    python3 same_arm_pairs.py reduce <file.root> -o <out_stem>
        -> <out_stem>.parquet  one row per (event, arm, lepton track) in a gap
        -> <out_stem>.json     per-file counts for the normalisation

Row columns:
  event, arm (0..3; 2 = A, 3 = C), track, parent, particle (e-/e+),
  op / ov        creator process / birth volume of the track
  par_op         creator process of the parent if the parent left a hit
  gamma_E        for a conversion lepton: the converting gamma's energy (MeV)
  ox, oy, oz     birth vertex (global mm)
  n_steps, edep_eV, w_min, w_max    its DriftGas steps (prompt)
  u0, v0, w0, gx0, gy0, gz0, px0, py0, pz0, ke0   first gap step (local / global
                 position, global momentum direction, kinetic energy MeV)
  wall_eV, plas_eV   its own prompt deposit in this arm's SiPM wall / plastic
  cap            capture category (thermal_accounting.capture_cat)
  ent            the primary entered the capsule envelope (gun r < 11.5 mm)
  trig_arms      bitmask of arms with a legacy trigger leg (prompt hits only)
  neutron_E_eV   primary energy

Choices: PROMPT = hit time < 1e8 ns (28Al beta decays are removed; their
gamma-less betas are kept only if prompt, i.e. never).  A lepton is "in the gap"
with >= 3 prompt DriftGas steps.  The trigger is thermal_accounting's legacy leg
(per arm: a SiPM bar >= 0.5 MIP and a plastic bar >= 0.5 MIP, no time window).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thermal_accounting as TA  # noqa: E402

BR = ["eventID", "trackID", "parentID", "armID", "detType", "particle", "u", "v", "w",
      "edep", "ke", "time", "gx", "gy", "gz", "px", "py", "pz", "origin_vol",
      "origin_proc", "ox", "oy", "oz"]
MIN_STEPS = 3


def reduce_file(fp: str, stem: str) -> None:
    import pandas as pd
    import uproot
    C = Counter()
    with uproot.open(fp) as f:
        E = TA.event_table(f)
        et = f["EventTree"].arrays(["eventID", "neutron_E_eV"], library="np")
        conv = (f["ConvPairTree"].arrays(["eventID", "gamma_E", "gamma_trackID", "conv_vol"],
                                         library="np") if "ConvPairTree" in f else None)
        parts = []
        for t in f["HitTree"].iterate(BR, library="np", step_size="200 MB"):
            det = t["detType"].astype(str)
            prt = t["particle"].astype(str)
            keep = (t["time"] < TA.T_PROMPT_NS) & (prt != "neutron") & (prt != "gamma") & \
                np.isin(det, ["DriftGas", "PlasticScint", "BackScintL", "BackScintR"])
            d = {k: t[k][keep] for k in BR if k not in ("detType", "particle", "origin_vol",
                                                         "origin_proc")}
            d["det"] = det[keep]
            d["prt"] = prt[keep]
            d["ov"] = np.char.rstrip(t["origin_vol"][keep].astype(str), "\x00")
            d["op"] = np.char.rstrip(t["origin_proc"][keep].astype(str), "\x00")
            parts.append(pd.DataFrame(d))
    H = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=BR)
    nE = np.zeros(E["emax"])
    nE[et["eventID"]] = et["neutron_E_eV"]
    C["N.sim"] = int(E["nev"])
    C["N.entering"] = int(E["ent"][E["eventID"]].sum())
    for k, name in enumerate(E["cat_names"]):
        C[f"cap.{name}"] = int(((E["cat"][E["eventID"]] == k) & E["ent"][E["eventID"]]).sum())

    # ---- trigger legs (legacy: per arm, a SiPM bar and a plastic bar >= 0.5 MIP)
    S = H[H.det != "DriftGas"].copy()
    S["chan"] = np.where(S.det == "PlasticScint",
                         np.clip((S.u + 250.0) // 25.0, 0, 19), np.where(S.det == "BackScintL", 20, 21))
    g = S.groupby(["eventID", "armID", "chan"]).edep.sum().reset_index()
    g["ok"] = np.where(g.chan < 20, g.edep >= TA.LEG_MIP * TA.MIP_S_EV, g.edep >= TA.LEG_MIP * TA.MIP_P_EV)
    g["is_p"] = g.chan >= 20
    a = g[g.ok].groupby(["eventID", "armID"]).is_p.agg(["any", "all"])
    legs = a[a["any"] & ~a["all"]].reset_index()[["eventID", "armID"]]
    trig = legs.groupby("eventID").armID.agg(lambda s: int(sum(1 << int(x) for x in s)))
    C["trig.events"] = int(len(trig))
    for arm in range(4):
        C[f"trig.legs.arm{arm}"] = int((legs.armID == arm).sum())
    C["trig.events.entering"] = int(E["ent"][trig.index.to_numpy()].sum()) if len(trig) else 0
    for k, name in enumerate(E["cat_names"]):
        ix = trig.index.to_numpy()
        C[f"trig.events.cap.{name}"] = int((E["cat"][ix] == k).sum()) if len(ix) else 0

    # ---- leptons in a drift gap
    D = H[(H.det == "DriftGas") & H.prt.isin(["e-", "e+"])]
    D = D.sort_values(["eventID", "armID", "trackID", "time"], kind="stable")
    key = ["eventID", "armID", "trackID"]
    agg = D.groupby(key).agg(n_steps=("edep", "size"), edep_eV=("edep", "sum"),
                             w_min=("w", "min"), w_max=("w", "max"))
    first = D.groupby(key).first()[["parentID", "prt", "op", "ov", "ox", "oy", "oz", "u", "v",
                                    "w", "gx", "gy", "gz", "px", "py", "pz", "ke"]]
    L = agg.join(first).reset_index()
    L = L[L.n_steps >= MIN_STEPS]
    # own scintillator deposit per (event, arm, track)
    sw = S[S.det == "PlasticScint"].groupby(key).edep.sum().rename("wall_eV")
    sp = S[S.det != "PlasticScint"].groupby(key).edep.sum().rename("plas_eV")
    L = L.join(sw, on=key).join(sp, on=key)
    L[["wall_eV", "plas_eV"]] = L[["wall_eV", "plas_eV"]].fillna(0.0)
    # parent's creator, where the parent left any hit in a scored volume
    info = H.groupby(["eventID", "trackID"]).op.first()
    pk = pd.MultiIndex.from_arrays([L.eventID.to_numpy(), L.parentID.to_numpy()])
    L["par_op"] = info.reindex(pk).to_numpy()
    # converting gamma
    L["gamma_E"] = np.nan
    if conv is not None and len(conv["eventID"]):
        cm = pd.Series(conv["gamma_E"], index=pd.MultiIndex.from_arrays(
            [conv["eventID"].astype(np.int64), conv["gamma_trackID"].astype(np.int64)]))
        cm = cm[~cm.index.duplicated()]
        isconv = L.op == "conv"
        L.loc[isconv, "gamma_E"] = cm.reindex(pd.MultiIndex.from_arrays(
            [L.eventID[isconv].astype(np.int64), L.parentID[isconv].astype(np.int64)])).to_numpy()
    ev = L.eventID.to_numpy()
    L["cap"] = np.asarray(E["cat_names"])[E["cat"][ev]]
    L["ent"] = E["ent"][ev]
    L["trig_arms"] = trig.reindex(ev).fillna(0).astype(int).to_numpy()
    L["neutron_E_eV"] = nE[ev]
    L = L.rename(columns={"eventID": "event", "armID": "arm", "trackID": "track",
                          "parentID": "parent", "prt": "particle", "u": "u0", "v": "v0",
                          "w": "w0", "gx": "gx0", "gy": "gy0", "gz": "gz0", "px": "px0",
                          "py": "py0", "pz": "pz0", "ke": "ke0"})
    L["file"] = os.path.basename(fp).replace("_t0.root", "")
    for c in L.columns:
        if L[c].dtype == np.float64 and c not in ("gamma_E",):
            L[c] = L[c].astype(np.float32)
    L.to_parquet(stem + ".parquet", index=False)
    C["rows"] = int(len(L))
    for arm in range(4):
        n = L[L.arm == arm].groupby("event").size()
        C[f"gap.events1.arm{arm}"] = int((n >= 1).sum())
        C[f"gap.events2.arm{arm}"] = int((n >= 2).sum())
    with open(stem + ".json", "w") as fo:
        json.dump(dict(schema="mx17/same_arm_pairs/1", file=fp, C=dict(C)), fo, indent=1)
    print(f"{fp}: {len(L)} gap leptons, {C['trig.events']} trigger events, "
          f"two-lepton arm-events A {C['gap.events2.arm2']} C {C['gap.events2.arm3']}")


def main() -> int:
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("reduce")
    r.add_argument("inp")
    r.add_argument("-o", "--out", required=True, help="output stem")
    a = ap.parse_args()
    reduce_file(a.inp, a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
