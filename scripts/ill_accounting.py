#!/usr/bin/env python3
"""
ill_accounting.py -- per-neutron background budget for the ILL cell scan
========================================================================
The §7 "background per absorbed neutron" rows of HANDOFF_SIM.md
(x17_facility_search/ill), from `mx17_full_sim --beam ill --target cell` runs.
Thresholds, MIP calibration and the trigger-leg definition are imported from
thermal_accounting.py so the ILL numbers and the n_TOF contract are comparable.

    # per file (lxplus, LCG_106) -> small JSON of weighted counts
    python3 ill_accounting.py reduce <file.root> -o part.json
    # merge -> accounting.json (+ one CSV row for the scan table)
    python3 ill_accounting.py merge parts/*.json -o <outdir> --config G1

Every count is weighted by EventTree.weight (1/bias for a biased ³He(n,γ),
≈1 otherwise), so it is per *unbiased* primary.  Divided by the (n,p)
absorptions in He3Gas it becomes "per absorbed neutron".

DEFINITIONS (written to the JSON):
* target volumes: He3Gas and every He3Cell_* volume (classified by logical
  volume, never by radius).
* gap charge: ≥ 1 keV deposited by charged tracks in one arm's DriftGas
  (thermal_accounting Tier A).  A reactor beam is continuous, so delayed
  deposits (activation β, e.g. ²⁸Al, T½ = 2.24 min) are a steady-state
  background at the same rate as prompt ones: both are reported, and "all"
  is the one to use.
* trigger leg: in one arm, a SiPM-wall bar ≥ 0.5 MIP and a plastic bar
  ≥ 0.5 MIP (charge summed per channel, the legacy event-level emulation of
  thermal_accounting F5); pair-tag: legs in ≥ 2 arms.  Prompt hits only
  (< 100 ms); a delayed leg pairs with nothing from the same neutron.
* wall pair in the gaps: a γ→e⁺e⁻ conversion (ConvPairTree) in a target
  volume whose e⁻ AND e⁺ (parentID = converting γ) both leave hits in
  DriftGas; "two arms" if they do so in different arms.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermal_accounting import (MIP_P_EV, MIP_S_EV, LEG_MIP, TIER_A_EV,  # noqa: E402
                                T_PROMPT_NS, NEUTRAL, SIPM_DET, PLAS_DET,
                                Vocab, per_group_sum, HIT_BRANCHES)
from collections import defaultdict  # noqa: E402


def read_hits(f, voc: Vocab):
    """thermal_accounting.read_hits plus global hit positions (gx, gy, gz)."""
    parts = defaultdict(list)
    if not f["HitTree"].num_entries:
        return None
    sid = None
    for t in f["HitTree"].iterate(HIT_BRANCHES + ["gx", "gy", "gz"], library="np",
                                  step_size="150 MB"):
        det = voc.enc("det", t["detType"])
        prt = voc.enc("prt", t["particle"])
        nid = voc.ids("prt", ["neutron"])
        sid = voc.ids("det", (SIPM_DET,) + PLAS_DET)
        keep = ~np.isin(prt, nid) | np.isin(det, sid)
        parts["det"].append(det[keep])
        parts["prt"].append(prt[keep])
        parts["ov"].append(voc.enc("ov", t["origin_vol"][keep]))
        parts["op"].append(voc.enc("op", t["origin_proc"][keep]))
        for k in ("eventID", "trackID", "parentID", "armID", "edep", "time", "u",
                  "ox", "oy", "oz", "gx", "gy", "gz"):
            parts[k].append(t[k][keep])
    h = {k: np.concatenate(v) for k, v in parts.items()}
    order = np.argsort(h["eventID"], kind="stable")
    return {k: v[order] for k, v in h.items()}

SCHEMA = "ill/accounting/1"
TARGET_PREFIX = ("He3Gas", "He3Cell_", "He3Cap_", "Slab")


def is_target(v: str) -> bool:
    return v.startswith(TARGET_PREFIX)


def s_(x) -> str:
    return (x.decode() if isinstance(x, bytes) else str(x)).rstrip("\x00")


def reduce_file(fp: str) -> dict:
    import uproot
    C = Counter()
    with uproot.open(fp) as f:
        et = f["EventTree"].arrays(["eventID", "capture_vol", "capture_proc", "weight",
                                    "n_thrown"], library="np")
        conv = f["ConvPairTree"].arrays(["eventID", "conv_vol", "gamma_trackID", "gamma_E"],
                                        library="np") if "ConvPairTree" in f else None
        voc = Vocab()
        h = read_hits(f, voc)

    ev = et["eventID"].astype(np.int64)
    w = et["weight"].astype(float)
    emax = int(ev.max()) + 1
    W = np.zeros(emax)
    W[ev] = w
    C["N.sim"] = int(len(ev))
    C["N.thrown"] = int(et["n_thrown"].sum())
    C["W.sum"] = float(w.sum())
    vols = np.array([s_(x) for x in et["capture_vol"]])
    procs = np.array([s_(x).replace("biasWrapper(", "").rstrip(")") for x in et["capture_proc"]])
    for (v, p) in set(zip(vols, procs)):
        m = (vols == v) & (procs == p)
        C[f"budget.{v or '(escaped)'}.{p or '-'}"] += float(w[m].sum())
        C[f"budget_raw.{v or '(escaped)'}.{p or '-'}"] += int(m.sum())
    C["absorbed_np"] = float(w[(vols == "He3Gas") & (procs == "neutronInelastic")].sum())

    if h is None:
        return dict(schema=SCHEMA, C=C, _table=None)

    V, names = voc.code, voc.names
    hev = h["eventID"].astype(np.int64)
    arm, det, prt = h["armID"], h["det"], h["prt"]
    edep, tm = h["edep"], h["time"]
    charged = ~np.isin(prt, voc.ids("prt", list(NEUTRAL)))
    prompt = tm < T_PROMPT_NS
    is_dg = det == V["det"].get("DriftGas", -9)
    is_s = det == V["det"].get(SIPM_DET, -9)
    is_p = np.isin(det, voc.ids("det", PLAS_DET))

    # ---- gap charge (Tier A) ---------------------------------------------
    def tier_a(mask):
        out = np.zeros(emax, bool)
        if mask.any():
            k, s = per_group_sum(hev[mask] * 4 + arm[mask], edep[mask])
            out[np.unique(k[s >= TIER_A_EV] // 4)] = True
        return out

    ta_p = tier_a(is_dg & charged & prompt)
    ta_all = tier_a(is_dg & charged)
    C["gap.prompt"] = float(W[ta_p].sum())
    C["gap.all"] = float(W[ta_all].sum())
    C["gap_raw.all"] = int(ta_all.sum())
    # where the gap charge comes from: origin volume of the largest deposit
    m = is_dg & charged
    if m.any():
        idx = np.flatnonzero(m)
        order = np.lexsort((-edep[idx], hev[idx]))
        idx = idx[order]
        first = np.r_[True, hev[idx][1:] != hev[idx][:-1]]
        OV, OP = names["ov"], names["op"]
        for i in idx[first]:
            e = hev[i]
            if not ta_all[e]:
                continue
            ov = OV[h["ov"][i]]
            cls = "target" if is_target(ov) else ("air" if ov == "World" else "detector/other")
            C[f"gap_origin.{cls}"] += float(W[e])
            C[f"gap_origin_proc.{OP[h['op'][i]]}"] += float(W[e])
            if cls == "target":
                C[f"gap_origin_target.{ov}"] += float(W[e])

    # ---- trigger legs / pair-tags (legacy, per channel) ------------------
    u = h["u"]
    chan = np.full(len(hev), -1, np.int64)
    chan[is_s] = np.clip((u[is_s] + 250.0) // 25.0, 0, 19).astype(np.int64)
    chan[det == V["det"].get("BackScintL", -9)] = 20
    chan[det == V["det"].get("BackScintR", -9)] = 21

    def legs(mask):
        if not mask.any():
            return np.array([], np.int64)
        k, s = per_group_sum((hev[mask] * 4 + arm[mask]) * 22 + chan[mask], edep[mask])
        s_ok = (k % 22 < 20) & (s >= LEG_MIP * MIP_S_EV)
        p_ok = (k % 22 >= 20) & (s >= LEG_MIP * MIP_P_EV)
        return np.intersect1d(np.unique(k[s_ok] // 22), np.unique(k[p_ok] // 22))

    for tag, mask in (("prompt", (is_s | is_p) & prompt), ("all", is_s | is_p)):
        ea = legs(mask)
        evl = ea // 4
        C[f"legs.{tag}"] = float(W[evl].sum())
        C[f"legs_raw.{tag}"] = int(len(ea))
        ue, cnt = np.unique(evl, return_counts=True)
        C[f"trig_events.{tag}"] = float(W[ue].sum())
        C[f"pairtags.{tag}"] = float(W[ue[cnt >= 2]].sum())
        C[f"pairtags_raw.{tag}"] = int((cnt >= 2).sum())

    # ---- per-event, per-arm table (any trigger menu / accidentals offline) --
    is_l = det == V["det"].get("LiqScint_1", -9)
    sc = is_s | is_p | is_l
    T = None
    if sc.any() or ta_all.any():
        evs = np.union1d(np.unique(hev[sc]), np.flatnonzero(ta_all))
        inn = np.isin(hev, evs)          # hits of events in the table only
        row = np.clip(np.searchsorted(evs, hev), 0, len(evs) - 1)
        T = dict(ev=evs, w=W[evs].astype(np.float32))
        z = lambda: np.zeros((len(evs), 4), np.float32)  # noqa: E731
        m = is_dg & charged & inn
        T["gap"] = z()
        np.add.at(T["gap"], (row[m], arm[m]), edep[m] * 1e-6)
        # edep-weighted gap centroid per arm, mm (global frame)
        T["gpos"] = np.full((len(evs), 4, 3), np.nan, np.float32)
        if m.any():
            g = np.zeros((len(evs), 4, 3))
            for c, key in enumerate(("gx", "gy", "gz")):
                np.add.at(g[:, :, c], (row[m], arm[m]), edep[m] * h[key][m])
            ok = T["gap"] > 0
            T["gpos"][ok] = (g[ok] / (T["gap"][ok] * 1e6)[:, None]).astype(np.float32)
        for tag, tm_ok in (("p", prompt), ("a", np.ones_like(prompt))):
            for key, mk in (("sipm", is_s), ("plast", is_p)):
                mm_ = mk & tm_ok
                T[f"{key}_{tag}"] = z()
                if mm_.any():
                    k, ss = per_group_sum((hev[mm_] * 4 + arm[mm_]) * 22 + chan[mm_], edep[mm_])
                    ea = k // 22
                    np.maximum.at(T[f"{key}_{tag}"], (np.searchsorted(evs, ea // 4), ea % 4),
                                  (ss * 1e-6).astype(np.float32))
            mm_ = is_l & tm_ok & inn
            T[f"ls_{tag}"] = z()
            np.add.at(T[f"ls_{tag}"], (row[mm_], arm[mm_]), edep[mm_] * 1e-6)
        T["t_scint"] = np.full((len(evs), 4), np.inf, np.float32)
        m = sc & prompt & inn
        np.minimum.at(T["t_scint"], (row[m], arm[m]), tm[m].astype(np.float32))
        cv = {v: i for i, v in enumerate(sorted(set(vols)))}
        T["capvol"] = np.array([cv[v] for v in vols[np.searchsorted(ev, evs)]], np.int16)
        T["capvol_names"] = np.array(sorted(cv, key=cv.get))

    # ---- wall pairs reaching the gaps ------------------------------------
    if conv is not None and len(conv["eventID"]):
        cvol = np.array([s_(x) for x in conv["conv_vol"]])
        ce = conv["eventID"].astype(np.int64)
        gt = conv["gamma_trackID"].astype(np.int64)
        tgt = np.array([is_target(v) for v in cvol])
        for v in set(cvol[tgt]):
            C[f"conv.in_target.{v}"] += float(W[ce[cvol == v]].sum())
        # DriftGas hits of conversion daughters: (event, parentID, particle, arm)
        md = is_dg & np.isin(prt, voc.ids("prt", ["e-", "e+"]))
        em_id = V["prt"].get("e-", -9)
        dg = {}
        for e, p, pt, a in zip(hev[md], h["parentID"][md], prt[md], arm[md]):
            dg.setdefault((int(e), int(p)), {}).setdefault(int(pt == em_id), set()).add(int(a))
        for e, g, v, ok in zip(ce, gt, cvol, tgt):
            if not ok:
                continue
            d = dg.get((int(e), int(g)))
            if d and 0 in d and 1 in d:
                C[f"wallpair_gaps.{v}"] += float(W[e])
                C["wallpair_gaps"] += float(W[e])
                C["wallpair_gaps_raw"] += 1
                if len(d[0] | d[1]) > 1:
                    C["wallpair_gaps_2arm"] += float(W[e])
    return dict(schema=SCHEMA, C=C, _table=T)


def merge(parts, outdir, config):
    C = Counter()
    n = 0
    for p in parts:
        d = json.load(open(p))
        C.update(d["C"])
        n += 1
    N = C["N.sim"]
    A = C["absorbed_np"]
    per_n = {k: v / N for k, v in C.items() if not k.startswith(("N.", "budget_raw", ))
             and not k.endswith("_raw") and "_raw." not in k}
    per_abs = {k: v / A for k, v in per_n.items() if A}
    per_abs = {k: v * N for k, v in per_abs.items()}
    out = dict(schema=SCHEMA, config=config, files=n, N_sim=N, N_thrown=C["N.thrown"],
               aperture_acceptance=N / C["N.thrown"] if C["N.thrown"] else None,
               absorbed_np_per_primary=A / N, counts=dict(C), per_primary=per_n,
               per_absorbed=per_abs)
    od = Path(outdir)
    od.mkdir(parents=True, exist_ok=True)
    (od / "accounting.json").write_text(json.dumps(out, indent=1, sort_keys=True))
    print(f"{config}: {n} parts, {N:,} primaries, absorbed (n,p) {A/N:.5f}/primary")
    for k in sorted(per_abs):
        if k.startswith(("budget.", "gap", "legs.", "pairtags.", "trig_events", "wallpair", "conv.")):
            print(f"  {k:55s} {per_abs[k]:.4e} /abs   raw {C.get(k.replace('.', '_raw.', 1), '')}")
    return out


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("reduce")
    r.add_argument("file")
    r.add_argument("-o", "--out", required=True)
    m = sp.add_parser("merge")
    m.add_argument("parts", nargs="+")
    m.add_argument("-o", "--outdir", required=True)
    m.add_argument("--config", default="")
    a = ap.parse_args()
    if a.cmd == "reduce":
        d = reduce_file(a.file)
        d["file"] = a.file
        T = d.pop("_table", None)
        if T is not None:
            np.savez_compressed(str(Path(a.out).with_suffix(".npz")), **T)
        Path(a.out).write_text(json.dumps(d))
    else:
        files = []
        for p in a.parts:
            files += sorted(glob.glob(p)) if "*" in p else [p]
        merge(files, a.outdir, a.config)


if __name__ == "__main__":
    main()
