#!/usr/bin/env python3
"""
thermal_accounting.py -- where the thermal neutrons go, step by step.

The funnel scoped in nTof_x17/ntof_athens_26/HANDOFF_THERMAL_ACCOUNTING.md:

    F1  what the neutron does            (per neutron entering the capsule)
    F2  does a capture put a charged particle in a Micromegas drift gap
    F3  where that charged particle was made (leading depositor)
    F4  which reaction made the capsule pair (+ internal pairs, analytic)
    F5  what fires a trigger leg (separate breakdown, all neutrons)

Two stages, so the heavy part runs where the data is:

    # per ROOT file (lxplus / condor, LCG_106) -> one small JSON of counts
    python3 thermal_accounting.py reduce  <file.root> -o part.json
    python3 thermal_accounting.py reduce  <file.root> -o part.json --mode pairs
    python3 thermal_accounting.py reduce  <file.root> -o part.json --mode bias

    # anywhere (numpy only) -> the contract: accounting.json + F1..F5.csv
    python3 thermal_accounting.py merge   parts/neutrons/*.json \\
        --pairs parts/pairs/*.json --bias parts/bias/*.json -o <outdir>

Inputs (nose-first campaign, 2026-07-23, MX17_Full_Geant 3d97437):
    neutrons_thermal_trig_2cm_nose     100 x 10M EAR2 neutrons, 1 meV - 2 eV
    pairs_thermal_trig_2cm_nose        X17 + IPC pairs at the gas capture vertices
    neutrons_thermal_bias1e5_2cm_nose  3He(n,g) biased x1e5 (cross-check only)

WHAT THE STORED TRUTH ALLOWS, AND THE CHOICES MADE (all written to the JSON):

* "Entering the capsule" is not stored.  The beam runs along +y from y = -200 mm
  and the capsule axis is the beam axis, so a primary with gun radius
  r = hypot(x, z) < 11.5 mm (the CFRP barrel) enters it on its straight line.
* Only the PRIMARY neutron's terminal interaction is stored (EventTree
  capture_vol / capture_proc).  In the thermal window that is the capture.
  ³He(n,p)t is "neutronInelastic" in He3Gas; radiative capture is "nCapture".
* Capture nucleus is stored only as a volume.  He3Cap_Al -> ²⁷Al.  He3Cap_CFRP
  is C + H + O (2.07 % H by mass in this geometry), so a CFRP capture is not
  split by nucleus here; the converting γ energy is histogrammed instead.
* Ancestry is one level per stored hit (trackID, parentID, creator process,
  birth volume and vertex).  A delta ray is walked up to its parent while the
  parent also left a hit; the ROOT of a deposit is the first track whose
  creator is not an ionisation process.  A chain that breaks (the parent left
  no hit in a scored volume) is kept as "unresolved" and counted.
* PROMPT window: hit time < 100 ms after the neutron is generated.  The sim
  also tracks radioactive decay (²⁸Al, T½ = 2.24 min) to completion, so about
  half of all drift-gas hits are β decays seconds to hours later.  In the data
  those are uncorrelated with the pulse; they are counted separately, never
  mixed into the prompt funnel.
* Tier A: >= 1 keV by charged tracks, summed per arm, in DriftGas.
  Tier B: >= 0.5 MIP by charged tracks in one SiPM-wall bar or plastic bar.
* Trigger leg (F5): in one arm, a SiPM-wall bar cluster >= 0.5 MIP and a plastic
  bar cluster >= 0.5 MIP within 20 ns (clusters: hits on one channel with gaps
  <= 20 ns; all particles, as the scintillator sees them).  MIP = 0.4751 /
  3.3494 MeV (2 cm muon MPV, 2026-07-26).  The legacy event-level emulation of
  analyze_trigger_thermal.py (no time window) is counted alongside.
* Leading depositor (F5): root carrying the most energy in the leg's two
  clusters; under 70 % of the leg's energy -> "mixed".
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np

SCHEMA = "athens/thermal_accounting/1"

CAPSULE = ("He3Gas", "He3Cap_Al", "He3Cap_CFRP")
R_ENVELOPE_MM = 11.5
R_BORE_MM = 10.0          # the gas bore: a straight path inside it crosses gas
R_NEAR_MM = 30.0
T_PROMPT_NS = 1.0e8
TIER_A_EV = 1.0e3
MIP_S_EV, MIP_P_EV = 0.4751e6, 3.3494e6
LEG_MIP = 0.5
TIER_B_MIP = 0.5
COINC_NS = 20.0
LEAD_SHARE = 0.70

CHAIN_PROCS = {"eIoni", "hIoni", "ionIoni", "muIoni"}
HADRONIC = {"hadElastic", "neutronInelastic", "nCapture", "protonInelastic",
            "dInelastic", "tInelastic", "alphaInelastic", "ionInelastic",
            "He3Inelastic", "neutronElastic"}
LEPTONS = {"e-", "e+"}
NEUTRAL = {"gamma", "neutron"}

SIPM_DET = "PlasticScint"
PLAS_DET = ("BackScintL", "BackScintR")

# normalisation
N_PULSE_SIM_WINDOW = 4.284e6     # EAR2 Ph3 flux integral, [1e-3, 2] eV, whole beam
PULSES_PER_DAY = 1.929e4         # x17_rate_3He.txt header
DAYS = 30
#: x17_rate_3He.txt, neutrons per pulse on the cell: 0.01-0.1 and 0.1-1 eV in
#: full, 1-10 eV up to 2 eV assuming iso-lethargic (log10 2 of the decade).  The
#: table has no bin below 10 meV.
RT_NEUTRONS_PER_PULSE = 3.11e6 + 1.68e6 + 1.39e6 * math.log10(2.0)

HIT_BRANCHES = ["eventID", "trackID", "parentID", "armID", "detType", "particle",
                "edep", "time", "u", "origin_vol", "origin_proc", "ox", "oy", "oz"]

CAP_CATS = ("he3_np", "ncapture_he3", "ncapture_wall_al", "ncapture_wall_cfrp",
            "capsule_other", "ncapture_elsewhere", "other_elsewhere",
            "no_absorption")
CAPTURE_CATS = ("ncapture_he3", "ncapture_wall_al", "ncapture_wall_cfrp",
                "ncapture_elsewhere")
NUCLEUS = {"ncapture_wall_al": "al27", "ncapture_wall_cfrp": "cfrp",
           "ncapture_he3": "he3", "ncapture_elsewhere": "capture_elsewhere"}


# --------------------------------------------------------------------------- #
# classification
# --------------------------------------------------------------------------- #
def s_(x) -> str:
    return x.decode() if isinstance(x, bytes) else str(x).rstrip("\x00")


def vol_class(v: str) -> str:
    if v in CAPSULE:
        return "capsule"
    if v in ("DriftGas", "AmpGas"):
        return "drift_gas"
    if v.startswith(("GasWindow", "DriftCathode", "Micromesh", "PCB_",
                     "ResistivePaste")):
        return "chamber"
    if v in ("PlasticScint", "BackScintL", "BackScintR", "LiqScint_1"):
        return "scintillator"
    if v.startswith(("BackScintAl", "BackScintTape", "LS_")):
        return "scint_wrap_vessel"
    if v == "World":
        return "air"
    if v == "":
        return "none"
    return "other"


def capture_cat(vol: str, proc: str) -> str:
    if proc.startswith("biasWrapper(") and proc.endswith(")"):
        proc = proc[len("biasWrapper("):-1]          # the biased campaign's name
    if vol == "":
        return "no_absorption"
    if vol == "He3Gas":
        if proc == "neutronInelastic":
            return "he3_np"
        if proc == "nCapture":
            return "ncapture_he3"
        return "capsule_other"
    if vol == "He3Cap_Al" and proc == "nCapture":
        return "ncapture_wall_al"
    if vol == "He3Cap_CFRP" and proc == "nCapture":
        return "ncapture_wall_cfrp"
    if vol in CAPSULE:
        return "capsule_other"
    if proc == "nCapture":
        return "ncapture_elsewhere"
    return "other_elsewhere"


def f3_bucket(op: str, ov: str, prt: str):
    """(bucket, sub) of a charged root track, by how and where it was made."""
    if op == "Radioactivation":
        return "other", "activation"
    if op in CHAIN_PROCS:
        return "other", "unresolved_ionisation"
    if prt not in LEPTONS:
        return "neutron_induced", prt if prt in ("proton", "deuteron", "alpha",
                                                   "triton") else "ion"
    if op == "conv":
        if ov in CAPSULE:
            return "pair_near_capsule", ov
        return "pair_elsewhere", vol_class(ov)
    if op in ("compt", "phot"):
        if ov in CAPSULE:
            return "compton_near_capsule", ov
        return "compton_elsewhere", vol_class(ov)
    if op == "nCapture":
        # conversion electrons emitted by G4PhotonEvaporation in the capture
        return "other", "capture_conversion_electron"
    if op in HADRONIC:
        return "neutron_induced", op
    return "other", op


def f5_bucket(root: dict, share: float, cat: str):
    """(bucket, sub) of a trigger leg from its leading root and the capture."""
    if share < LEAD_SHARE:
        return "mixed", ""
    op, ov, prt = root["op"], root["ov"], root["prt"]
    if op == "Radioactivation":
        return "activation", vol_class(ov)
    b, sub = f3_bucket(op, ov, prt) if prt != "gamma" else ("gamma_direct", "")
    if b == "pair_near_capsule":
        return b, NUCLEUS.get(cat, cat)
    if b == "neutron_induced":
        return b, sub
    if cat not in ("ncapture_wall_al", "ncapture_wall_cfrp", "ncapture_he3"):
        return "capture_outside_capsule", cat
    if b in ("compton_near_capsule", "compton_elsewhere", "gamma_direct"):
        return "compton_capsule_capture", NUCLEUS[cat]
    if b == "pair_elsewhere":
        return "conv_structure_capsule_capture", sub
    return "other", sub


# --------------------------------------------------------------------------- #
# reading
# --------------------------------------------------------------------------- #
class Vocab:
    """String columns -> small ints, consistent across chunks of one file."""

    def __init__(self):
        self.code = defaultdict(dict)
        self.names = defaultdict(list)

    def enc(self, field: str, arr) -> np.ndarray:
        u, inv = np.unique(np.asarray(arr).astype(str), return_inverse=True)
        d, n = self.code[field], self.names[field]
        lut = np.empty(len(u), np.int32)
        for i, s in enumerate(u):
            s = s.rstrip("\x00")
            if s not in d:
                d[s] = len(n)
                n.append(s)
            lut[i] = d[s]
        return lut[inv]

    def ids(self, field: str, names) -> np.ndarray:
        d = self.code[field]
        return np.array([d[x] for x in names if x in d], np.int32)


def event_table(f):
    et = f["EventTree"].arrays(["eventID", "capture_vol", "capture_proc",
                                "vtx_x", "vtx_z", "weight"], library="np")
    nev = len(et["eventID"])
    emax = int(et["eventID"].max()) + 1
    cv_u, cv_i = np.unique(et["capture_vol"].astype(str), return_inverse=True)
    cp_u, cp_i = np.unique(et["capture_proc"].astype(str), return_inverse=True)
    combo = cv_i * len(cp_u) + cp_i
    cu, ci = np.unique(combo, return_inverse=True)
    cat_names = list(CAP_CATS)
    sub_names = ["capsule", "drift_gas", "chamber", "scintillator",
                 "scint_wrap_vessel", "air", "none", "other"]
    lut_cat = np.empty(len(cu), np.int8)
    lut_sub = np.empty(len(cu), np.int8)
    for k, c in enumerate(cu):
        vol = cv_u[c // len(cp_u)].rstrip("\x00")
        proc = cp_u[c % len(cp_u)].rstrip("\x00")
        lut_cat[k] = cat_names.index(capture_cat(vol, proc))
        lut_sub[k] = sub_names.index(vol_class(vol))
    ev_cat = np.full(emax, -1, np.int8)
    ev_sub = np.full(emax, -1, np.int8)
    ev_ent = np.zeros(emax, bool)
    ev_bore = np.zeros(emax, bool)
    r = np.hypot(et["vtx_x"], et["vtx_z"])
    ev_cat[et["eventID"]] = lut_cat[ci]
    ev_sub[et["eventID"]] = lut_sub[ci]
    ev_ent[et["eventID"]] = r < R_ENVELOPE_MM
    ev_bore[et["eventID"]] = r < R_BORE_MM
    return dict(nev=nev, emax=emax, cat=ev_cat, sub=ev_sub, ent=ev_ent, bore=ev_bore,
                cat_names=cat_names, sub_names=sub_names,
                eventID=et["eventID"], weight=et["weight"])


def read_hits(f, voc: Vocab, drop_neutrons=True):
    keep_det = (SIPM_DET,) + PLAS_DET
    parts = defaultdict(list)
    if not f["HitTree"].num_entries:
        return None
    for t in f["HitTree"].iterate(HIT_BRANCHES, library="np", step_size="150 MB"):
        det = voc.enc("det", t["detType"])
        prt = voc.enc("prt", t["particle"])
        keep = np.ones(len(det), bool)
        if drop_neutrons:
            nid = voc.ids("prt", ["neutron"])
            sid = voc.ids("det", keep_det)
            keep = ~np.isin(prt, nid) | np.isin(det, sid)
        parts["det"].append(det[keep])
        parts["prt"].append(prt[keep])
        parts["ov"].append(voc.enc("ov", t["origin_vol"][keep]))
        parts["op"].append(voc.enc("op", t["origin_proc"][keep]))
        for k in ("eventID", "trackID", "parentID", "armID", "edep", "time", "u",
                  "ox", "oy", "oz"):
            parts[k].append(t[k][keep])
    h = {k: np.concatenate(v) for k, v in parts.items()}
    order = np.argsort(h["eventID"], kind="stable")
    return {k: v[order] for k, v in h.items()}


# --------------------------------------------------------------------------- #
# reduce: neutrons
# --------------------------------------------------------------------------- #
def per_group_sum(keys: np.ndarray, w: np.ndarray):
    u, inv = np.unique(keys, return_inverse=True)
    return u, np.bincount(inv, weights=w, minlength=len(u))


def reduce_neutrons(fp: str) -> dict:
    import uproot
    C = Counter()
    Hh = defaultdict(Counter)
    voc = Vocab()
    with uproot.open(fp) as f:
        E = event_table(f)
        conv = f["ConvPairTree"].arrays(["eventID", "conv_vol", "gamma_E",
                                         "gamma_trackID"], library="np") \
            if "ConvPairTree" in f else None
        h = read_hits(f, voc)

    cat, sub, ent = E["cat"], E["sub"], E["ent"]
    CN, SN = E["cat_names"], E["sub_names"]
    ev_all = E["eventID"]
    C["meta.nonunit_weight"] = int((E["weight"] != 1.0).sum())
    C["N.sim"] = E["nev"]
    C["N.entering"] = int(ent[ev_all].sum())
    C["N.bore"] = int(E["bore"][ev_all].sum())

    # ---- F1 ------------------------------------------------------------- #
    for k, name in enumerate(CN):
        m = cat[ev_all] == k
        C[f"F1.{name}"] = int((m & ent[ev_all]).sum())
        C[f"F1bore.{name}"] = int((m & E["bore"][ev_all]).sum())
        C[f"miss.{name}"] = int((m & ~ent[ev_all]).sum())
        if name in ("ncapture_elsewhere", "other_elsewhere"):
            for j, sname in enumerate(SN):
                n = int((m & ent[ev_all] & (sub[ev_all] == j)).sum())
                if n:
                    C[f"F1.{name}.{sname}"] = n
                n = int((m & ~ent[ev_all] & (sub[ev_all] == j)).sum())
                if n:
                    C[f"miss.{name}.{sname}"] = n

    # ---- conversions produced in the capsule (for the internal-pair acceptance)
    gmap = {}
    if conv is not None and len(conv["eventID"]):
        cvol = conv["conv_vol"].astype(str)
        incap = np.isin(np.char.rstrip(cvol, "\x00"), CAPSULE)
        ce = conv["eventID"]
        for e, g, tid, v in zip(ce, conv["gamma_E"], conv["gamma_trackID"], cvol):
            gmap[(int(e), int(tid))] = (float(g), v.rstrip("\x00"))
        evc = np.unique(ce[incap])
        for k, name in enumerate(CN):
            sel = evc[(cat[evc] == k) & ent[evc]]
            if len(sel):
                C[f"conv.capsule_events.{name}"] += int(len(sel))
        for e, g in zip(ce[incap], conv["gamma_E"][incap]):
            if ent[e]:
                Hh[f"conv.capsule_gammaE.{CN[cat[e]]}"][int(g / 0.25)] += 1
        C["conv.capsule_n.entering"] = int(ent[ce[incap]].sum())

    if h is None:
        return dict(C=C, H=Hh)

    V = voc.code
    names = voc.names
    ev, arm, det, prt = h["eventID"], h["armID"], h["det"], h["prt"]
    edep, tm, u = h["edep"], h["time"], h["u"]
    charged = ~np.isin(prt, voc.ids("prt", list(NEUTRAL)))
    prompt = tm < T_PROMPT_NS
    is_dg = det == V["det"].get("DriftGas", -9)
    is_s = det == V["det"].get(SIPM_DET, -9)
    is_p = np.isin(det, voc.ids("det", PLAS_DET))
    chan = np.full(len(ev), -1, np.int64)
    chan[is_s] = np.clip((u[is_s] + 250.0) // 25.0, 0, 19).astype(np.int64)
    chan[det == V["det"].get("BackScintL", -9)] = 20
    chan[det == V["det"].get("BackScintR", -9)] = 21
    ev64 = ev.astype(np.int64)
    emax = E["emax"]

    # time structure of drift-gas charged hits, for the record
    radio = h["op"] == V["op"].get("Radioactivation", -9)
    lt = np.log10(np.maximum(tm[is_dg & charged], 1e-3))
    for x, r in zip(np.floor(lt * 2).astype(int), radio[is_dg & charged]):
        Hh["time.dg_charged_log10ns_x2." + ("radioactivation" if r else "other")][int(x)] += 1

    def arm_flag(mask, thr):
        out = np.zeros(emax, bool)
        if mask.any():
            k, s = per_group_sum(ev64[mask] * 4 + arm[mask], edep[mask])
            out[np.unique(k[s >= thr] // 4)] = True
        return out

    tierA = arm_flag(is_dg & charged & prompt, TIER_A_EV)
    tierA_del = arm_flag(is_dg & charged & ~prompt, TIER_A_EV)
    tierB = np.zeros(emax, bool)
    mb = (is_s | is_p) & charged & prompt
    if mb.any():
        k, s = per_group_sum((ev64[mb] * 4 + arm[mb]) * 22 + chan[mb], edep[mb])
        thr = np.where(k % 22 < 20, TIER_B_MIP * MIP_S_EV, TIER_B_MIP * MIP_P_EV)
        tierB[np.unique(k[s >= thr] // 88)] = True

    # ---- F2 ------------------------------------------------------------- #
    for name in CN:
        k = CN.index(name)
        sel = ev_all[(cat[ev_all] == k) & ent[ev_all]]
        pre = "F2" if name in CAPTURE_CATS else "F2x"
        C[f"{pre}.{name}.total"] = int(len(sel))
        C[f"{pre}.{name}.tierA"] = int(tierA[sel].sum())
        C[f"{pre}.{name}.tierB"] = int(tierB[sel].sum())
        C[f"{pre}.{name}.tierA_or_B"] = int((tierA[sel] | tierB[sel]).sum())
        C[f"{pre}.{name}.tierA_delayed_only"] = int((tierA_del[sel] & ~tierA[sel]).sum())
        sel_m = ev_all[(cat[ev_all] == k) & ~ent[ev_all]]
        C[f"F2miss.{name}.total"] = int(len(sel_m))
        C[f"F2miss.{name}.tierA"] = int(tierA[sel_m].sum())

    # ---- per-event hit slices ------------------------------------------- #
    ue, starts = np.unique(ev64, return_index=True)
    ends = np.r_[starts[1:], len(ev64)]
    slice_of = dict(zip(ue.tolist(), zip(starts.tolist(), ends.tolist())))

    OP, OV, PRT = names["op"], names["ov"], names["prt"]
    trk_l, par_l = h["trackID"], h["parentID"]

    def ancestry(a, b):
        info = {}
        tr, pa = trk_l[a:b].tolist(), par_l[a:b].tolist()
        ops, ovs, prs = h["op"][a:b].tolist(), h["ov"][a:b].tolist(), prt[a:b].tolist()
        ox, oy, oz = h["ox"][a:b].tolist(), h["oy"][a:b].tolist(), h["oz"][a:b].tolist()
        for i in range(len(tr)):
            if tr[i] not in info:
                info[tr[i]] = dict(par=pa[i], op=OP[ops[i]], ov=OV[ovs[i]],
                                   prt=PRT[prs[i]], r=math.sqrt(ox[i] ** 2 + oy[i] ** 2 + oz[i] ** 2),
                                   trk=tr[i])
        memo = {}

        def root(t):
            if t in memo:
                return memo[t]
            seen, cur = [], t
            while True:
                seen.append(cur)
                d = info[cur]
                if d["op"] in CHAIN_PROCS and d["par"] in info and len(seen) < 64:
                    cur = d["par"]
                    continue
                break
            for s in seen:
                memo[s] = cur
            return cur
        return info, root

    # ---- F3 / F4 -------------------------------------------------------- #
    capk = [CN.index(c) for c in CAPTURE_CATS]
    cand = ev_all[np.isin(cat[ev_all], capk) & ent[ev_all] & tierA[ev_all]]
    m_dg = is_dg & charged & prompt
    for e in cand.tolist():
        a, b = slice_of[e]
        info, root = ancestry(a, b)
        Eroot = defaultdict(float)
        for i in np.flatnonzero(m_dg[a:b]).tolist():
            Eroot[root(int(trk_l[a + i]))] += float(edep[a + i])
        if not Eroot:
            continue
        lead = max(Eroot, key=Eroot.get)
        share = Eroot[lead] / sum(Eroot.values())
        d = info[lead]
        bucket, sb = f3_bucket(d["op"], d["ov"], d["prt"])
        cname = CN[cat[e]]
        C["F3.total"] += 1
        C[f"F3.{bucket}"] += 1
        C[f"F3.{bucket}.{sb}"] += 1
        C[f"F3.by_capture.{cname}.{bucket}"] += 1
        Hh["F3.lead_share_x10"][min(int(share * 10), 9)] += 1
        if d["op"] in ("conv", "compt", "phot") and d["r"] < R_NEAR_MM:
            C[f"F3.sens_r30.{'pair' if d['op'] == 'conv' else 'compton'}_near_capsule"] += 1
        if bucket == "pair_near_capsule":
            nuc = NUCLEUS[cname]
            C[f"F4.external.{nuc}"] += 1
            g = gmap.get((e, d["par"]))
            if g is not None:
                Hh[f"F4.gammaE_x4.{nuc}"][int(g[0] * 4)] += 1
                if g[0] <= 2.3:
                    C[f"F4.external.{nuc}.gammaE_le_2p3MeV"] += 1
            else:
                C[f"F4.external.{nuc}.gamma_not_found"] += 1

    # ---- F5 ------------------------------------------------------------- #
    msp = is_s | is_p
    if msp.any():
        k, s = per_group_sum((ev64[msp] * 4 + arm[msp]) * 22 + chan[msp], edep[msp])
        s_ok = (k % 22 < 20) & (s >= LEG_MIP * MIP_S_EV)
        p_ok = (k % 22 >= 20) & (s >= LEG_MIP * MIP_P_EV)
        ea_s, ea_p = np.unique(k[s_ok] // 22), np.unique(k[p_ok] // 22)
        ea = np.intersect1d(ea_s, ea_p)             # (event, arm) legacy legs
        C["F5.legacy.legs"] = int(len(ea))
        ev_leg = ea // 4
        C["F5.legacy.trigger_events"] = int(len(np.unique(ev_leg)))
        uu, cc = np.unique(ev_leg, return_counts=True)
        C["F5.legacy.pairtags"] = int((cc >= 2).sum())
    else:
        ea = np.array([], np.int64)

    by_event = defaultdict(list)
    for x in ea.tolist():
        by_event[x // 4].append(x % 4)
    for e, arms in by_event.items():
        a, b = slice_of[e]
        info, root = ancestry(a, b)
        legs = []
        for am in arms:
            idx = a + np.flatnonzero(msp[a:b] & (arm[a:b] == am))
            clus = []                              # (is_plastic, E, t0, [idx])
            for ch in np.unique(chan[idx]).tolist():
                ii = idx[chan[idx] == ch]
                ii = ii[np.argsort(tm[ii], kind="stable")]
                t_ = tm[ii]
                cut = np.flatnonzero(np.diff(t_) > COINC_NS) + 1
                for grp in np.split(ii, cut):
                    clus.append((ch >= 20, float(edep[grp].sum()), float(tm[grp[0]]), grp))
            best = None
            for cs in clus:
                if cs[0] or cs[1] < LEG_MIP * MIP_S_EV:
                    continue
                for cp in clus:
                    if not cp[0] or cp[1] < LEG_MIP * MIP_P_EV:
                        continue
                    if abs(cs[2] - cp[2]) > COINC_NS:
                        continue
                    if best is None or (cp[1], cs[1]) > (best[1][1], best[0][1]):
                        best = (cs, cp)
            if best is None:
                continue
            t_leg = min(best[0][2], best[1][2])
            hits = np.r_[best[0][3], best[1][3]]
            Er = defaultdict(float)
            for i in hits.tolist():
                Er[root(int(trk_l[i]))] += float(edep[i])
            lead = max(Er, key=Er.get)
            share = Er[lead] / sum(Er.values())
            bucket, sb = f5_bucket(info[lead], share, CN[cat[e]])
            pd_ = "prompt" if t_leg < T_PROMPT_NS else "delayed"
            C[f"F5.legs.{pd_}"] += 1
            C[f"F5.legs.{pd_}.{bucket}"] += 1
            if sb:
                C[f"F5.legs.{pd_}.{bucket}.{sb}"] += 1
            if ent[e]:
                C[f"F5.legs_entering.{pd_}"] += 1
            legs.append(dict(arm=am, t=t_leg, lead=lead, info=info[lead],
                             bucket=bucket, sub=sb, Ep=best[1][1], prompt=pd_ == "prompt"))
        pl = sorted([l for l in legs if l["prompt"]], key=lambda l: -l["Ep"])
        if len({l["arm"] for l in pl}) >= 2:
            l1 = pl[0]
            l2 = next(l for l in pl[1:] if l["arm"] != l1["arm"])
            C["F5.pairtags.prompt"] += 1
            dt20 = abs(l1["t"] - l2["t"]) <= COINC_NS
            if dt20:
                C["F5.pairtags.prompt_dt20"] += 1
            if l1["lead"] == l2["lead"]:
                rel = "one_particle"
            elif (l1["info"]["op"] == "conv" and l2["info"]["op"] == "conv"
                  and l1["info"]["par"] == l2["info"]["par"]):
                rel = "genuine_pair"
            else:
                rel = "two_ancestries_same_capture"
            bk = l1["bucket"] if l1["bucket"] == l2["bucket"] else "mixed_buckets"
            C[f"F5.pairtags.prompt.relation.{rel}"] += 1
            C[f"F5.pairtags.prompt.bucket.{bk}"] += 1
            if dt20:
                C[f"F5.pairtags.prompt_dt20.relation.{rel}"] += 1
                C[f"F5.pairtags.prompt_dt20.bucket.{bk}"] += 1
    return dict(C=C, H=Hh)


# --------------------------------------------------------------------------- #
# reduce: pairs (acceptance of 20.6 MeV IPC/X17 pairs) and bias (3He(n,g))
# --------------------------------------------------------------------------- #
def reduce_pairs(fp: str) -> dict:
    import uproot
    C = Counter()
    with uproot.open(fp) as f:
        et = f["EventTree"].arrays(["eventID", "event_type"], library="np")
        h = f["HitTree"].arrays(["eventID", "trackID", "armID", "detType",
                                 "particle", "edep"], library="np")
    typ = {0: "X17", 1: "IPC"}
    emax = int(et["eventID"].max()) + 1
    ety = np.full(emax, -9)
    ety[et["eventID"]] = et["event_type"]
    for k, n in typ.items():
        C[f"pairs.{n}.n_gen"] = int((et["event_type"] == k).sum())
    det = np.char.rstrip(h["detType"].astype(str), "\x00")
    prt = np.char.rstrip(h["particle"].astype(str), "\x00")
    m = (det == "DriftGas") & ~np.isin(prt, list(NEUTRAL))
    ev = h["eventID"][m].astype(np.int64)
    k, s = per_group_sum(ev * 4 + h["armID"][m], h["edep"][m])
    evA = np.unique(k[s >= TIER_A_EV] // 4)
    # both primary leptons (trackID 1, 2) leave >= 1 keV in some drift gap
    mp = m & np.isin(h["trackID"], (1, 2))
    k2, s2 = per_group_sum(h["eventID"][mp].astype(np.int64) * 3 + h["trackID"][mp],
                           h["edep"][mp])
    ok = k2[s2 >= TIER_A_EV]
    ue, cnt = np.unique(ok // 3, return_counts=True)
    evB = ue[cnt >= 2]
    for kk, n in typ.items():
        C[f"pairs.{n}.tierA"] = int((ety[evA] == kk).sum())
        C[f"pairs.{n}.both_leptons_tierA"] = int((ety[evB] == kk).sum())
    return dict(C=C, H={})


def reduce_bias(fp: str) -> dict:
    import uproot
    C = Counter()
    with uproot.open(fp) as f:
        et = f["EventTree"].arrays(["capture_vol", "capture_proc", "vtx_x",
                                    "vtx_z", "weight"], library="np")
    vol = np.char.rstrip(et["capture_vol"].astype(str), "\x00")
    proc = np.char.rstrip(et["capture_proc"].astype(str), "\x00")
    proc = np.where(proc == "biasWrapper(nCapture)", "nCapture", proc)
    ent = np.hypot(et["vtx_x"], et["vtx_z"]) < R_ENVELOPE_MM
    m = (vol == "He3Gas") & (proc == "nCapture")
    C["bias.N.sim"] = int(len(vol))
    C["bias.N.entering"] = int(ent.sum())
    C["bias.he3_ng.mc_count"] = int((m & ent).sum())
    C["bias.he3_ng.sum_w"] = float(et["weight"][m & ent].sum())
    C["bias.he3_ng.sum_w2"] = float((et["weight"][m & ent] ** 2).sum())
    C["bias.he3_np.entering"] = int(((vol == "He3Gas") & (proc == "neutronInelastic") & ent).sum())
    return dict(C=C, H={})


# --------------------------------------------------------------------------- #
# merge -> the contract
# --------------------------------------------------------------------------- #
#: Analytic inputs, not computed here.  Values from nTof_x17/sept26_prelim_analysis
#: (ipc_aluminium, ipc_born), evaluated 2026-09-15; see provenance in the JSON.
ANALYTIC = {
    "pairs_per_wall_capture_M1": 2.1035e-3,     # capsule_summary('M1'), Al+C weighted
    "pairs_per_wall_capture_E1": 2.6654e-3,     # capsule_summary('E1'), the bracket
    "pairs_per_al_capture_M1": 2.097e-3,
    "pairs_per_c_capture_M1": 2.158e-3,
    "pairs_per_h_capture_M1": 4.516e-4,         # ipc_born.alpha_pair('M1', 2.2246 MeV), one γ per capture
    "pairs_per_he3_radiative_capture": 4.6646e-3,   # ipc_born M1 + E0
    "he3_ng_per_neutron": 1.0311e-8,            # bookkeeping(0.031), self-shielded
    "source": "nTof_x17 sept26_prelim_analysis.ipc_aluminium (capsule_summary, "
              "bookkeeping(0.031 eV), rate_comparison attrs) and ipc_born; NOT "
              "IPC/capture = 2.1e-3 of Viviani et al. (MeV).",
}


def load_parts(paths):
    C, H = Counter(), defaultdict(Counter)
    for p in paths:
        d = json.load(open(p))
        for k, v in d["C"].items():
            C[k] += v
        for name, hist in d.get("H", {}).items():
            for b, n in hist.items():
                H[name][int(b)] += n
    return C, H, len(paths)


def binom_err(n, N):
    if not N:
        return None
    p = n / N
    return math.sqrt(max(p * (1 - p), 0.0) / N)


def merge(args) -> int:
    parts = sorted(sum((glob.glob(p) for p in args.parts), []))
    C, H, nfiles = load_parts(parts)
    CP, _, npf = load_parts(sorted(sum((glob.glob(p) for p in args.pairs or []), [])))
    CB, _, nbf = load_parts(sorted(sum((glob.glob(p) for p in args.bias or []), [])))

    N_sim, N_ent = C["N.sim"], C["N.entering"]
    f_ent = N_ent / N_sim
    rt_factor = RT_NEUTRONS_PER_PULSE * PULSES_PER_DAY * DAYS / N_ent       # per MC entering neutron
    sim_factor = N_PULSE_SIM_WINDOW * PULSES_PER_DAY * DAYS / N_sim        # per MC neutron
    nodes = []

    def node(id_, parent, label, mc, parent_mc, source="geant4", note="",
             per_neutron=None, slide=True):
        pn = per_neutron if per_neutron is not None else (mc / N_ent if mc is not None else None)
        frac = (mc / parent_mc) if (mc is not None and parent_mc) else None
        nodes.append(dict(
            id=id_, parent=parent, label=label, mc_count=mc, weight=mc,
            frac_of_parent=frac,
            frac_err=binom_err(mc, parent_mc) if (mc is not None and parent_mc) else None,
            per_neutron=pn,
            per30d_ratetable=pn * RT_NEUTRONS_PER_PULSE * PULSES_PER_DAY * DAYS if pn is not None else None,
            per30d_sim=(mc * sim_factor if mc is not None else pn * f_ent * N_PULSE_SIM_WINDOW * PULSES_PER_DAY * DAYS),
            source=source, note=note, on_slide=slide))

    # ---------------- F1 ----------------
    node("root", None, "neutrons entering the capsule (E_n < 2 eV)", N_ent, None,
         note=f"gun radius < {R_ENVELOPE_MM} mm; {f_ent:.4f} of simulated neutrons")
    f1 = {c: C[f"F1.{c}"] for c in CAP_CATS}
    ng_all = f1["ncapture_wall_al"] + f1["ncapture_wall_cfrp"] + f1["ncapture_he3"] + f1["ncapture_elsewhere"]
    other = f1["capsule_other"] + f1["other_elsewhere"]
    node("F1.he3_np", "root", "³He(n,p) → p + t", f1["he3_np"], N_ent,
         note="p 573 keV and t 191 keV stop in the gas: invisible")
    node("F1.ncapture", "root", "(n,γ) capture, anywhere", ng_all, N_ent)
    node("F1.other", "root", "other", other, N_ent,
         note="non-capture terminal interactions outside the gas, e.g. ¹⁴N(n,p) in air")
    node("F1.no_absorption", "root", "leaves the world", f1["no_absorption"], N_ent)
    node("F1.ncapture.wall_al", "F1.ncapture", "²⁷Al(n,γ), capsule vessel",
         f1["ncapture_wall_al"], ng_all)
    node("F1.ncapture.wall_cfrp", "F1.ncapture", "CFRP wrap (n,γ): ¹H / ¹²C / ¹⁶O",
         f1["ncapture_wall_cfrp"], ng_all,
         note="nucleus not stored; geometry CFRP is 89.7 % C, 2.07 % H, 8.3 % O by mass")
    node("F1.ncapture.he3", "F1.ncapture", "³He(n,γ)⁴He", f1["ncapture_he3"], ng_all,
         source="analytic", per_neutron=ANALYTIC["he3_ng_per_neutron"],
         note=f"MC count {f1['ncapture_he3']} not quoted; analytic self-shielded at 31 meV "
              f"(ipc_aluminium.bookkeeping)")
    node("F1.ncapture.elsewhere", "F1.ncapture", "captured outside the capsule after scattering",
         f1["ncapture_elsewhere"], ng_all)
    for k, v in sorted(C.items()):
        if k.startswith("F1.ncapture_elsewhere."):
            node(f"F1.ncapture.elsewhere.{k.split('.')[-1]}", "F1.ncapture.elsewhere",
                 k.split(".")[-1], v, f1["ncapture_elsewhere"], slide=False)

    # ---------------- F2 ----------------
    tot = sum(C[f"F2.{c}.total"] for c in CAPTURE_CATS)
    chA = sum(C[f"F2.{c}.tierA"] for c in CAPTURE_CATS)
    chB = sum(C[f"F2.{c}.tierB"] for c in CAPTURE_CATS)
    delA = sum(C[f"F2.{c}.tierA_delayed_only"] for c in CAPTURE_CATS)
    node("F2.gamma_only", "F1.ncapture", "γ only: nothing charged in a drift gap (prompt)",
         tot - chA, tot)
    node("F2.charged_in_detector", "F1.ncapture",
         "≥ 1 charged particle in a Micromegas drift gap (≥ 1 keV, prompt)", chA, tot)
    node("F2.tierB.charged_in_scint", "F1.ncapture",
         "tier B: ≥ 0.5 MIP charged in a SiPM-wall or plastic bar (prompt)", chB, tot, slide=False)
    node("F2.delayed_only", "F1.ncapture",
         "drift-gap charge only from delayed decays (≥ 100 ms, ²⁸Al β)", delA, tot, slide=False,
         note="not in the prompt funnel; uncorrelated with the pulse in the data")
    f1_id = dict(ncapture_wall_al="wall_al", ncapture_wall_cfrp="wall_cfrp",
                 ncapture_he3="he3", ncapture_elsewhere="elsewhere")
    for c in CAPTURE_CATS:
        n, a_ = C[f"F2.{c}.total"], C[f"F2.{c}.tierA"]
        node(f"F2.by.{NUCLEUS[c]}.charged", f"F1.ncapture.{f1_id[c]}",
             f"{NUCLEUS[c]}: charged in drift gap", a_, n, slide=False)

    # ---------------- F3 ----------------
    f3tot = C["F3.total"]
    for b, lab in (("pair_near_capsule", "e⁺e⁻ pair made in the capsule"),
                   ("compton_near_capsule", "Compton / photo-electron made in the capsule"),
                   ("pair_elsewhere", "pair made elsewhere"),
                   ("compton_elsewhere", "Compton / photo-electron made elsewhere"),
                   ("neutron_induced", "neutron-induced (recoils, (n,p), capture recoils)"),
                   ("other", "other")):
        node(f"F3.{b}", "F2.charged_in_detector", lab, C[f"F3.{b}"], f3tot)
        for k, v in sorted(C.items()):
            if k.startswith(f"F3.{b}.") and not k.startswith("F3.by_"):
                node(f"F3.{b}.{k.split('.', 2)[2]}", f"F3.{b}", k.split(".", 2)[2], v,
                     C[f"F3.{b}"], slide=False)
    if C["F3.total"] != chA:
        nodes[-1]["note"] += f" [F3.total {C['F3.total']} vs tier A {chA}]"

    # ---------------- F4 ----------------
    # internal pairs: analytic per capture x the detector acceptance the sim
    # measures for EXTERNAL capsule pairs from the same capture volume.
    acc = {}
    for c in ("ncapture_wall_al", "ncapture_wall_cfrp"):
        made = C[f"conv.capsule_events.{c}"]
        acc[c] = C[f"F4.external.{NUCLEUS[c]}"] / made if made else None
    ipc_acc = CP["pairs.IPC.tierA"] / CP["pairs.IPC.n_gen"] if CP.get("pairs.IPC.n_gen") else None
    int_al = f1["ncapture_wall_al"] / N_ent * ANALYTIC["pairs_per_al_capture_M1"] * (acc["ncapture_wall_al"] or 0)
    # CFRP: this geometry's binder makes ~96 % of CFRP captures ¹H(n,γ); its single
    # 2.22 MeV M1 line gives 4.5e-4 pairs, which outweighs carbon's 3.7 % share
    # at 2.2e-3.  ¹⁶O is negligible.  The acceptance is the Al one: CFRP external
    # pairs are too few to measure their own, and the CFRP is the Al's outer skin.
    sh_c, sh_h = C_SHARE_OF_CFRP_CAPTURES, H_SHARE_OF_CFRP_CAPTURES
    int_cf_c = f1["ncapture_wall_cfrp"] / N_ent * sh_c * ANALYTIC["pairs_per_c_capture_M1"] * (acc["ncapture_wall_al"] or 0)
    int_cf_h = f1["ncapture_wall_cfrp"] / N_ent * sh_h * ANALYTIC["pairs_per_h_capture_M1"] * (acc["ncapture_wall_al"] or 0)
    int_cf = int_cf_c + int_cf_h
    int_he = ANALYTIC["he3_ng_per_neutron"] * ANALYTIC["pairs_per_he3_radiative_capture"] * (ipc_acc or 0)
    ext = {n: C[f"F4.external.{n}"] for n in ("al27", "cfrp", "he3", "capture_elsewhere")}
    ext_pn = sum(ext.values()) / N_ent
    parent_pn = ext_pn + int_al + int_cf + int_he
    node("F4.capsule_pairs", "F3.pair_near_capsule",
         "capsule e⁺e⁻ pairs with ≥ 1 lepton in a drift gap: external (Geant4) + internal (analytic)",
         None, None, source="geant4+analytic", per_neutron=parent_pn)

    def f4(id_, lab, pn, mc=None, source="geant4", note=""):
        nodes.append(dict(id=id_, parent="F4.capsule_pairs", label=lab, mc_count=mc, weight=mc,
                          frac_of_parent=pn / parent_pn if parent_pn else None,
                          frac_err=(math.sqrt(mc) / N_ent / parent_pn if (mc and parent_pn) else None),
                          per_neutron=pn,
                          per30d_ratetable=pn * RT_NEUTRONS_PER_PULSE * PULSES_PER_DAY * DAYS,
                          per30d_sim=pn * f_ent * N_PULSE_SIM_WINDOW * PULSES_PER_DAY * DAYS,
                          source=source, note=note, on_slide=True))

    f4("F4.external.al27", "²⁷Al(n,γ) → γ converts (external)", ext["al27"] / N_ent, ext["al27"])
    f4("F4.external.cfrp", "CFRP (n,γ) → γ converts (external)", ext["cfrp"] / N_ent, ext["cfrp"],
       note=f"{C['F4.external.cfrp.gammaE_le_2p3MeV']} of {ext['cfrp']} from a γ ≤ 2.3 MeV (¹H-like)")
    f4("F4.external.he3", "³He(n,γ) → γ converts (external)", ext["he3"] / N_ent, ext["he3"],
       note="MC; ≈ 0 expected")
    f4("F4.external.capture_elsewhere", "γ from a capture elsewhere, converting in the capsule",
       ext["capture_elsewhere"] / N_ent, ext["capture_elsewhere"])
    f4("F4.internal.al27", "²⁷Al(n,γ) internal pair creation", int_al, source="analytic",
       note=f"sim wall-Al captures x {ANALYTIC['pairs_per_al_capture_M1']:.3g}/capture (M1; all-E1 "
            f"x{ANALYTIC['pairs_per_wall_capture_E1'] / ANALYTIC['pairs_per_wall_capture_M1']:.2f}) "
            f"x acceptance {acc['ncapture_wall_al']} borrowed from external Al pairs")
    f4("F4.internal.cfrp", "CFRP (n,γ) internal pair creation: ¹H + ¹²C", int_cf, source="analytic",
       note=f"sim CFRP captures x (H share {sh_h:.3f} x {ANALYTIC['pairs_per_h_capture_M1']:.3g} "
            f"[{int_cf_h / int_cf if int_cf else 0:.0%}] + C share {sh_c:.3f} x "
            f"{ANALYTIC['pairs_per_c_capture_M1']:.3g}) x the Al acceptance {acc['ncapture_wall_al']}")
    f4("F4.internal.he3", "³He(n,γ) internal pair creation (IPC)", int_he, source="analytic",
       note=f"{ANALYTIC['he3_ng_per_neutron']:.3g}/n x {ANALYTIC['pairs_per_he3_radiative_capture']:.3g} "
            f"x IPC drift-gap acceptance {ipc_acc} (pairs_thermal_trig_2cm_nose, {npf} files)")

    # ---------------- F5 ----------------
    per_pulse = N_PULSE_SIM_WINDOW / N_sim
    legs = C["F5.legs.prompt"]
    f5nodes = []

    def f5(id_, parent, lab, mc, parent_mc, note="", slide=True):
        f5nodes.append(dict(id=id_, parent=parent, label=lab, mc_count=mc, weight=mc,
                            frac_of_parent=mc / parent_mc if parent_mc else None,
                            frac_err=binom_err(mc, parent_mc),
                            per_pulse=mc * per_pulse, per_neutron=mc / N_sim,
                            per30d_ratetable=mc / N_ent * RT_NEUTRONS_PER_PULSE * PULSES_PER_DAY * DAYS,
                            per30d_sim=mc * sim_factor, source="geant4", note=note, on_slide=slide))

    f5("F5.legs", None, "trigger legs, prompt (wall ∧ plastic, 0.5 MIP, 20 ns)", legs, None,
       note=f"legacy event-level emulation: {C['F5.legacy.legs'] * per_pulse:.1f} legs/pulse "
            f"(no time window); delayed legs {C['F5.legs.delayed'] * per_pulse:.2f}/pulse")
    f5("F5.legs_delayed", None, "trigger legs from delayed decays (≥ 100 ms)", C["F5.legs.delayed"], None,
       slide=False)
    for k, v in sorted(C.items()):
        p = k.split(".")
        if p[:3] == ["F5", "legs", "prompt"] and len(p) == 4:
            f5(f"F5.{p[3]}", "F5.legs", p[3], v, legs)
        elif p[:3] == ["F5", "legs", "prompt"] and len(p) == 5 and p[4]:
            f5(f"F5.{p[3]}.{p[4]}", f"F5.{p[3]}", p[4], v, C[f"F5.legs.prompt.{p[3]}"], slide=False)
    pt = C["F5.pairtags.prompt"]
    f5("F5.pairtags", None, "pair-tags: prompt legs in two arms, one event", pt, None,
       note=f"with |Δt| ≤ 20 ns: {C['F5.pairtags.prompt_dt20'] * per_pulse:.3f}/pulse. One neutron per "
            "event: the sim cannot make an accidental from two neutrons.")
    for k, v in sorted(C.items()):
        p = k.split(".")
        if p[:4] == ["F5", "pairtags", "prompt", "relation"]:
            f5(f"F5.pairtags.relation.{p[4]}", "F5.pairtags", p[4], v, pt)
        if p[:4] == ["F5", "pairtags", "prompt", "bucket"]:
            f5(f"F5.pairtags.bucket.{p[4]}", "F5.pairtags", p[4], v, pt, slide=False)

    out = dict(
        schema=SCHEMA,
        provenance=dict(
            git=args.git, campaign="neutrons_thermal_trig_2cm_nose", n_files=nfiles,
            n_primaries=N_sim, physics_list="FTFP_BERT_HP-style modular: EmStandard_option4, "
            "HadronElasticPhysicsHP, HadronPhysicsFTFP_BERT_HP, Decay, RadioactiveDecay",
            geant4=args.geant4, geometry="final surveyed 4-arm, nose-first capsule (3d97437)",
            window="E_n in [1 meV, 2 eV] (the >1 ms gate at 19.5 m)",
            near_capsule=dict(logical_volumes=list(CAPSULE), sensitivity=f"r < {R_NEAR_MM} mm from the capsule origin"),
            entering=f"gun radius < {R_ENVELOPE_MM} mm (beam along +y onto the capsule axis)",
            prompt_window_ns=T_PROMPT_NS, tierA_threshold_keV=TIER_A_EV / 1e3,
            tierB_threshold_MIP=TIER_B_MIP,
            trigger=dict(sipm_bar_MIP=LEG_MIP, plastic_bar_MIP=LEG_MIP, MIP_MeV=[MIP_S_EV / 1e6, MIP_P_EV / 1e6],
                         coincidence_ns=COINC_NS, cluster_gap_ns=COINC_NS, lead_share=LEAD_SHARE),
            ipc_added_analytically=True, internal_pairs_in_geant4=False,
            analytic=ANALYTIC, pairs_files=npf, bias_files=nbf,
            gaps=["capture nucleus inside CFRP not stored (C/H/O)",
                  "ancestry through tracks that left no scored hit is not followed",
                  "neutrons that enter, scatter out and are captured elsewhere are in F1.ncapture.elsewhere; "
                  "neutrons that miss the envelope but scatter in are not in the denominator"]),
        normalisation=dict(
            neutrons_simulated=N_sim, neutrons_entering_capsule=N_ent, entering_fraction=f_ent,
            pulses_per_day_ratetable=PULSES_PER_DAY,
            ratetable_neutrons_per_pulse_on_cell=RT_NEUTRONS_PER_PULSE,
            sim_neutrons_per_pulse_window=N_PULSE_SIM_WINDOW,
            per30d_ratetable_factor=rt_factor, per30d_sim_factor=sim_factor,
            per_pulse_sim_factor=per_pulse),
        nodes=nodes, trigger_nodes=f5nodes,
        cross_checks=dict(
            he3_ng_per_neutron_bias1e5=(CB["bias.he3_ng.sum_w"] / CB["bias.N.entering"]
                                        if CB.get("bias.N.entering") else None),
            he3_ng_bias1e5_mc_count=CB.get("bias.he3_ng.mc_count"),
            he3_ng_per_neutron_bias1e5_err=(math.sqrt(CB["bias.he3_ng.sum_w2"]) / CB["bias.N.entering"]
                                            if CB.get("bias.N.entering") else None),
            pairs=dict(CP), acceptance_external_capsule_pairs=acc,
            f3_sensitivity_r30=dict(pair=C["F3.sens_r30.pair_near_capsule"],
                                    compton=C["F3.sens_r30.compton_near_capsule"]),
            legacy_trigger=dict(legs_per_pulse=C["F5.legacy.legs"] * per_pulse,
                                trigger_events_per_pulse=C["F5.legacy.trigger_events"] * per_pulse,
                                pairtags_per_pulse=C["F5.legacy.pairtags"] * per_pulse),
            f1_gas_bore=dict(neutrons=C["N.bore"],
                             **{c: (C[f"F1bore.{c}"] / C["N.bore"] if C["N.bore"] else None)
                                for c in CAP_CATS},
                             note=f"gun radius < {R_BORE_MM} mm: straight path through the gas; "
                                  "compare ipc_aluminium.bookkeeping's 99.97 % (n,p)"),
            f2_missed_capsule={k: v for k, v in C.items() if k.startswith("F2miss.")},
            f2_non_capture={k: v for k, v in C.items() if k.startswith("F2x.")}),
        histograms={k: dict(sorted(v.items())) for k, v in H.items()},
        counts=dict(sorted(C.items())),
    )
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "accounting.json"), "w") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    cols = ["id", "parent", "label", "mc_count", "weight", "frac_of_parent", "frac_err",
            "per_neutron", "per_pulse", "per30d_ratetable", "per30d_sim", "source", "note", "on_slide"]
    groups = defaultdict(list)
    for n in nodes:
        groups[n["id"].split(".")[0] if n["id"] != "root" else "F1"].append(n)
    groups["F5"] = f5nodes
    for g, rows in groups.items():
        with open(os.path.join(args.out, f"{g}.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
    print("wrote", os.path.join(args.out, "accounting.json"), "and", ", ".join(f"{g}.csv" for g in groups))
    return 0


#: Share of CFRP captures on carbon, from this geometry's CFRP (mass fractions
#: C 0.8968, H 0.0207, O 0.0826) and thermal capture cross sections
#: σ(H) = 0.3326 b, σ(C) = 3.53 mb, σ(O) = 0.19 mb.
def _cfrp_shares():
    nC, nH, nO = 0.8968 / 12.011, 0.0207 / 1.008, 0.0826 / 15.999
    r = dict(C=nC * 3.53e-3, H=nH * 0.3326, O=nO * 1.9e-4)
    return r["C"] / sum(r.values()), r["H"] / sum(r.values())


C_SHARE_OF_CFRP_CAPTURES, H_SHARE_OF_CFRP_CAPTURES = _cfrp_shares()


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description="thermal-neutron accounting funnel")
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("reduce")
    r.add_argument("file")
    r.add_argument("-o", "--out", required=True)
    r.add_argument("--mode", choices=("neutrons", "pairs", "bias"), default="neutrons")
    m = sp.add_parser("merge")
    m.add_argument("parts", nargs="+")
    m.add_argument("--pairs", nargs="*")
    m.add_argument("--bias", nargs="*")
    m.add_argument("-o", "--out", required=True)
    m.add_argument("--git", default="3d97437")
    m.add_argument("--geant4", default="11.2 (LCG, lxplus build)")
    a = ap.parse_args()
    if a.cmd == "merge":
        return merge(a)
    t0 = time.time()
    fn = dict(neutrons=reduce_neutrons, pairs=reduce_pairs, bias=reduce_bias)[a.mode]
    res = fn(a.file)
    out = dict(schema=SCHEMA + "/part", mode=a.mode, file=a.file,
               elapsed_s=round(time.time() - t0, 1),
               C=dict(res["C"]), H={k: dict(v) for k, v in res["H"].items()})
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"{a.mode}: {a.file} -> {a.out}  ({out['elapsed_s']} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
