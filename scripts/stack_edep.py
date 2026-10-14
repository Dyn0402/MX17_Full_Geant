#!/usr/bin/env python3
"""
stack_edep.py -- per (event, arm) with a lepton in the drift gap: the TOTAL prompt
energy in each SiPM-wall strip and each plastic bar of that arm, for the
nTof_x17 scintillator-calorimetry check (ntof_calorimetry/wall_plastic.py).

same_arm_pairs.py keeps only each lepton's OWN deposit (its trackID), which
misses its delta rays, its bremsstrahlung, and the other particles of the event
in the same counter -- the hardware sees all of them.  This keeps the sums.

    python3 stack_edep.py reduce <file.root> -o <out_stem>
        -> <out_stem>.parquet   one row per (event, arm) with >= 1 gap lepton
        -> <out_stem>.json      per-file counts

Row columns:
  event, arm (0..3; 2 = A, 3 = C)
  w00..w19        wall ("PlasticScint", 3 mm) per 25 mm strip in u, MeV
  pL, pR          plastic bars (BackScintL / R), MeV
  n_lep           leptons in this arm's gap (>= 3 prompt steps)
  lead_*          the lepton with the most gap steps: n_steps, w_span (mm),
                  u0, v0, w0 (first gap step, local), px0, py0, pz0 (global
                  direction), ke0 (MeV), particle, op, ov,
                  own_wall, own_plas (its own deposits, MeV)
  cap, ent, neutron_E_eV
Prompt = hit time < thermal_accounting.T_PROMPT_NS (100 ms).  The hardware integrates
~160 ns around the trigger, so w*/pL/pR are summed only over hits within
WIN_NS of the leading lepton's first gap step; the 100 ms totals are kept as
wall_all / pL_all / pR_all for comparison.
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

BR = ["eventID", "trackID", "armID", "detType", "particle", "u", "v", "w",
      "edep", "ke", "time", "px", "py", "pz", "origin_vol", "origin_proc"]
MIN_STEPS = 3
WIN_NS = (-100.0, 200.0)


def reduce_file(fp: str, stem: str) -> None:
    import pandas as pd
    import uproot
    C = Counter()
    with uproot.open(fp) as f:
        E = TA.event_table(f)
        et = f["EventTree"].arrays(["eventID", "neutron_E_eV"], library="np")
        parts = []
        for t in f["HitTree"].iterate(BR, library="np", step_size="200 MB"):
            det = t["detType"].astype(str)
            prt = t["particle"].astype(str)
            keep = (t["time"] < TA.T_PROMPT_NS) & \
                np.isin(det, ["DriftGas", "PlasticScint", "BackScintL", "BackScintR"])
            d = {k: t[k][keep] for k in BR if k not in ("detType", "particle", "origin_vol",
                                                         "origin_proc")}
            d["det"] = det[keep]
            d["prt"] = prt[keep]
            d["ov"] = np.char.rstrip(t["origin_vol"][keep].astype(str), "\x00")
            d["op"] = np.char.rstrip(t["origin_proc"][keep].astype(str), "\x00")
            parts.append(pd.DataFrame(d))
    H = pd.concat(parts, ignore_index=True)
    nE = np.zeros(E["emax"])
    nE[et["eventID"]] = et["neutron_E_eV"]
    C["N.sim"] = int(E["nev"])

    # ---- leptons in a gap; the leading one per (event, arm)
    D = H[(H.det == "DriftGas") & H.prt.isin(["e-", "e+"])]
    D = D.sort_values(["eventID", "armID", "trackID", "time"], kind="stable")
    key = ["eventID", "armID", "trackID"]
    agg = D.groupby(key).agg(n_steps=("edep", "size"), w_min=("w", "min"), w_max=("w", "max"))
    first = D.groupby(key).first()[["prt", "op", "ov", "u", "v", "w", "px", "py", "pz", "ke", "time"]]
    L = agg.join(first).reset_index()
    L = L[L.n_steps >= MIN_STEPS]
    if not len(L):
        pd.DataFrame().to_parquet(stem + ".parquet")
        json.dump(dict(schema="mx17/stack_edep/2", file=fp, C=dict(C)), open(stem + ".json", "w"))
        return
    L["w_span"] = L.w_max - L.w_min
    n_lep = L.groupby(["eventID", "armID"]).size().rename("n_lep")
    lead = L.sort_values(["eventID", "armID", "n_steps"], ascending=[True, True, False]) \
        .groupby(["eventID", "armID"]).first()

    # ---- scintillator sums (every particle) and the leading lepton's own share
    S = H[H.det != "DriftGas"].copy()
    t_lead = lead["time"].rename("t_lead")
    S = S.join(t_lead, on=["eventID", "armID"])
    S_all = S.copy()
    S = S[(S.time - S.t_lead >= WIN_NS[0]) & (S.time - S.t_lead <= WIN_NS[1])]
    S_all["isw"] = S_all.det == "PlasticScint"
    allsum = S_all.groupby(["eventID", "armID", "det"]).edep.sum().unstack("det", fill_value=0.0) / 1e6
    allsum = allsum.reindex(columns=["PlasticScint", "BackScintL", "BackScintR"], fill_value=0.0)
    allsum.columns = ["wall_all", "pL_all", "pR_all"]
    S["chan"] = np.where(S.det == "PlasticScint", np.clip((S.u + 250.0) // 25.0, 0, 19),
                         np.where(S.det == "BackScintL", 20, 21)).astype(int)
    tot = S.groupby(["eventID", "armID", "chan"]).edep.sum().unstack("chan", fill_value=0.0)
    tot = tot.reindex(columns=range(22), fill_value=0.0) / 1e6
    tot.columns = [f"w{c:02d}" for c in range(20)] + ["pL", "pR"]
    own = S.groupby(["eventID", "armID", "trackID", "det"]).edep.sum()

    R = lead.join(n_lep).join(tot, how="left").join(allsum, how="left")
    R[list(tot.columns) + list(allsum.columns)] = R[list(tot.columns) + list(allsum.columns)].fillna(0.0)
    idx = R.index
    tk = R.trackID.to_numpy()
    ow = own.reindex(pd.MultiIndex.from_arrays([idx.get_level_values(0), idx.get_level_values(1), tk,
                                                np.full(len(R), "PlasticScint")])).fillna(0).to_numpy()
    opL = own.reindex(pd.MultiIndex.from_arrays([idx.get_level_values(0), idx.get_level_values(1), tk,
                                                 np.full(len(R), "BackScintL")])).fillna(0).to_numpy()
    opR = own.reindex(pd.MultiIndex.from_arrays([idx.get_level_values(0), idx.get_level_values(1), tk,
                                                 np.full(len(R), "BackScintR")])).fillna(0).to_numpy()
    R["lead_own_wall"] = ow / 1e6
    R["lead_own_plas"] = (opL + opR) / 1e6
    R = R.rename(columns={"n_steps": "lead_n_steps", "w_span": "lead_w_span", "prt": "lead_particle",
                          "op": "lead_op", "ov": "lead_ov", "u": "lead_u0", "v": "lead_v0",
                          "w": "lead_w0", "px": "lead_px0", "py": "lead_py0", "pz": "lead_pz0",
                          "ke": "lead_ke0", "time": "lead_t_ns"}).drop(columns=["trackID", "w_min", "w_max"])
    R = R.reset_index().rename(columns={"eventID": "event", "armID": "arm"})
    ev = R.event.to_numpy()
    R["cap"] = np.asarray(E["cat_names"])[E["cat"][ev]]
    R["ent"] = E["ent"][ev]
    R["neutron_E_eV"] = nE[ev]
    R["file"] = os.path.basename(fp).replace("_t0.root", "")
    for c in R.columns:
        if R[c].dtype == np.float64:
            R[c] = R[c].astype(np.float32)
    R.to_parquet(stem + ".parquet", index=False)
    C["rows"] = int(len(R))
    json.dump(dict(schema="mx17/stack_edep/2", file=fp, C=dict(C)), open(stem + ".json", "w"), indent=1)
    print(f"{fp}: {len(R)} (event, arm) rows")


def main() -> int:
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("reduce")
    r.add_argument("file")
    r.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    if a.cmd == "reduce":
        reduce_file(a.file, a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
