#!/usr/bin/env python3
"""digest_plastic_opt.py — per-(event,arm) trigger observables for optimizing the
PLASTIC threshold to trigger on real MM tracks (IPC) and reject double-Compton.

Per surviving (in-gate) event & arm:
  S   = max SiPM-bar edep            [MeV]   (SiPM-wall discriminator)
  P   = max plastic-bar (L/R) edep   [MeV]   (plastic discriminator = the knob)
  MM  = DriftGas+AmpGas edep by charged particles [MeV] — the TRACK label
        (a real e± leaves a track in the TPC gas; a double-Compton γ does not)
  EP  = e+ edep in this arm's SiPM+plastic [MeV] (>0 ⇒ pair-production, not Compton)
  Ssum= sum SiPM edep, Psum = sum plastic edep (diagnostics)

Signal (pairs) rows carry event_type (0=X17,1=IPC); background (neutrons) rows
carry isAl (capture_vol startswith He3Cap_Al). Thermal gate: TOF+t_hit>1 ms.

Run on lxplus. Writes npz to /afs/cern.ch/work/d/dneff/convpair_out/.
    python3 scripts/digest_plastic_opt.py signal|thermal [NFILES]
"""
import sys, glob
import numpy as np, uproot

WHICH = sys.argv[1] if len(sys.argv) > 1 else "signal"
NF = int(sys.argv[2]) if len(sys.argv) > 2 else 20
DIRS = {"signal": "/eos/experiment/ntof/data/x17/full_sim/pairs_thermal_trig_2cm",
        "thermal": "/eos/experiment/ntof/data/x17/full_sim/neutrons_thermal_trig_2cm"}
OUT = f"/eos/experiment/ntof/data/x17/full_sim/convpair_out/plastic_opt_{WHICH}.npz"
FLIGHT_MS, GATE_MS = 1.41, 1.0

def code_of(det):
    # det: object array of bytes (Char_t[32]/C). Raw == avoids costly decode.
    c = np.full(len(det), -1, np.int8)
    c[(det == b"PlasticScint") | (det == "PlasticScint")] = 0   # SiPM bar
    c[(det == b"BackScintL") | (det == "BackScintL")] = 1        # plastic L
    c[(det == b"BackScintR") | (det == "BackScintR")] = 2        # plastic R
    c[(det == b"DriftGas") | (det == "DriftGas") |
      (det == b"AmpGas") | (det == "AmpGas")] = 4                # MM
    return c

def digest(fp, is_signal):
    f = uproot.open(fp)
    evb = ["eventID", "event_type", "neutron_E_eV"]
    e = f["EventTree"].arrays(evb + (["capture_vol"] if not is_signal else []),
                              library="np")
    nev = f["EventTree"].num_entries
    mx = int(e["eventID"].max()) + 1
    eN = np.zeros(mx); eN[e["eventID"]] = e["neutron_E_eV"]
    typ = np.full(mx, -9, np.int64); typ[e["eventID"]] = e["event_type"]
    isAl = np.zeros(mx, bool)
    if not is_signal:
        cv = np.array([x.decode() if isinstance(x, bytes) else x for x in e["capture_vol"]]).astype(str)
        isAl[e["eventID"]] = np.char.startswith(cv, "He3Cap_Al")
    # accumulate per (event,arm) key = ev*4+arm; vectorized per-chunk reduction
    from collections import defaultdict
    barmax_s = defaultdict(float); barmax_p = defaultdict(float)
    mm = defaultdict(float); ep = defaultdict(float)
    ssum = defaultdict(float); psum = defaultdict(float)

    def merge_max(dst, keys, vals):
        for k, v in zip(keys.tolist(), vals.tolist()):
            if v > dst[k]: dst[k] = v
    def merge_add(dst, keys, vals):
        for k, v in zip(keys.tolist(), vals.tolist()):
            dst[k] += v

    for t in f["HitTree"].iterate(
            ["eventID", "armID", "detType", "edep", "time", "u", "particle"],
            library="np", step_size="250 MB"):
        c = code_of(t["detType"])
        keep = c >= 0
        e_hit = eN[t["eventID"]]
        with np.errstate(divide="ignore"):
            tof = np.where(e_hit > 0, FLIGHT_MS / np.sqrt(e_hit), np.inf)
        keep &= (tof + t["time"] * 1e-6) > GATE_MS
        if not keep.any():
            continue
        c = c[keep]; ev = t["eventID"][keep].astype(np.int64)
        arm = t["armID"][keep].astype(np.int64); ed = t["edep"][keep] * 1e-6
        u = t["u"][keep]
        praw = t["particle"][keep]     # object array of bytes; compare raw
        is_gamma = (praw == b"gamma") | (praw == "gamma")
        is_neu   = (praw == b"neutron") | (praw == "neutron")
        is_ep    = (praw == b"e+") | (praw == "e+")
        akey = ev * 4 + arm
        chan = np.where(c == 0, np.clip((u + 250.) // 25., 0, 19).astype(np.int64),
                        np.where(c == 1, 20, np.where(c == 2, 21, 30)))
        # bar maxima: per (ev,arm,chan) sum, then reduce to per-(ev,arm) max
        bkey = akey * 40 + chan
        ub, inv = np.unique(bkey, return_inverse=True)
        bar = np.zeros(len(ub)); np.add.at(bar, inv, ed)
        bch = ub % 40; bao = ub // 40
        sm = bch <= 19
        if sm.any():
            uk, uinv = np.unique(bao[sm], return_inverse=True)
            v = np.zeros(len(uk)); np.maximum.at(v, uinv, bar[sm]); merge_max(barmax_s, uk, v)
        pm = (bch == 20) | (bch == 21)
        if pm.any():
            uk, uinv = np.unique(bao[pm], return_inverse=True)
            v = np.zeros(len(uk)); np.maximum.at(v, uinv, bar[pm]); merge_max(barmax_p, uk, v)
        # summed observables: MM(charged), EP(e+ in scint), Ssum, Psum
        charged = ~(is_gamma | is_neu)
        for mask, dst in ((( c == 4) & charged, mm),
                          ((((c == 0)|(c == 1)|(c == 2)) & is_ep), ep),
                          ((c == 0), ssum),
                          (((c == 1)|(c == 2)), psum)):
            if mask.any():
                uk, uinv = np.unique(akey[mask], return_inverse=True)
                v = np.zeros(len(uk)); np.add.at(v, uinv, ed[mask]); merge_add(dst, uk, v)
    keys = set(barmax_s) | set(barmax_p)
    rows = []
    for key in keys:
        evid = key // 4
        rows.append((evid, key % 4, barmax_s.get(key, 0.), barmax_p.get(key, 0.),
                     mm.get(key, 0.), ep.get(key, 0.), ssum.get(key, 0.),
                     psum.get(key, 0.), typ[evid], isAl[evid]))
    return rows, nev

def main():
    files = sorted(glob.glob(f"{DIRS[WHICH]}/*.root"))[:NF]
    is_signal = WHICH == "signal"
    allrows = []; tot_ev = 0
    for k, fp in enumerate(files):
        try:
            r, nev = digest(fp, is_signal)
        except Exception as ex:
            print("  skip", fp.split("/")[-1], repr(ex)[:50]); continue
        allrows += r; tot_ev += nev
        print(f"  [{k+1}/{len(files)}] {fp.split('/')[-1]} rows={len(r)} cum={len(allrows)}", flush=True)
    A = np.array(allrows, dtype=object)
    ev = A[:, 0].astype(np.int64); arm = A[:, 1].astype(np.int64)
    S = A[:, 2].astype(float); P = A[:, 3].astype(float); MM = A[:, 4].astype(float)
    EP = A[:, 5].astype(float); Ss = A[:, 6].astype(float); Ps = A[:, 7].astype(float)
    typ = A[:, 8].astype(np.int64); isAl = A[:, 9].astype(bool)
    np.savez(OUT, ev=ev, arm=arm, S=S, P=P, MM=MM, EP=EP, Ssum=Ss, Psum=Ps,
             typ=typ, isAl=isAl, n_events=tot_ev)
    print(f"\n{WHICH}: events={tot_ev:,} arm-rows={len(ev):,}  wrote {OUT}")

if __name__ == "__main__":
    main()
