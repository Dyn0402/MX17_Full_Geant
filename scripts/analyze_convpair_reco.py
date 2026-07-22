#!/usr/bin/env python3
"""analyze_convpair_reco.py — Task 3 + 4c: detector-level (multiple-scattering)
opening angle of the Al conv-pair background, reconstructed from the MM DriftGas
hits, straight out of Geant4 (no analytic MSC model).

For each Al γ→e+e- conversion (ConvPairTree, conv_vol==He3Cap_Al), the daughter
leptons' MM hits are `HitTree` DriftGas hits with parentID == gamma_trackID,
split e-/e+ by the `particle` field. Each lepton's MM-ENTRY direction is the
truth momentum direction (px,py,pz) of its earliest DriftGas hit — i.e. the
direction AFTER scattering through the Al capsule + air, which is what a tracker
sees. The reconstructed opening angle is the angle between the two entry dirs.

Definitions:
  - "MM single"   : lepton has >=1 DriftGas hit  (reached the tracker)
  - "MM pair"     : BOTH leptons have DriftGas hits (reconstructable 2-track pair)
  - same-arm pair : both leptons enter the same MM arm (hard 2-track separation)

Rates (Task 3) use weight/pulse = FLUX_WINDOW / N_sim, ×PULSES_DAY.

Run on lxplus.  Optional arg: a single .root file (for validation), else the
whole neutrons_convpair_2cm campaign (completed files only).
"""
import sys, glob, os, json
import numpy as np, uproot

FLUX_WINDOW = 4_292_400.81
PULSES_DAY  = 1.929e4
DIR = "/eos/experiment/ntof/data/x17/full_sim/neutrons_convpair_2cm"
OUT = os.environ.get("CONVPAIR_RECO_OUT",
                     "/afs/cern.ch/work/d/dneff/convpair_out/convpair_reco.npz")

def first_dir(times, px, py, pz):
    """MM-entry unit direction = momentum dir of earliest-time hit."""
    j = int(np.argmin(times))
    v = np.array([px[j], py[j], pz[j]]); n = np.linalg.norm(v)
    return v / n if n > 0 else None

def process(fp, rows):
    f = uproot.open(fp)
    if "ConvPairTree" not in [k.split(";")[0] for k in f.keys()]:
        return 0
    c = f["ConvPairTree"].arrays(library="np")
    cv = c["conv_vol"]  # object array of bytes; compare without decoding
    al = np.where((cv == b"He3Cap_Al") | (cv == "He3Cap_Al"))[0]
    if len(al) == 0:
        return 0
    # event -> {gamma_trackID -> pair-index}
    ev_g = {}
    for i in al:
        ev_g.setdefault(int(c["eventID"][i]), {})[int(c["gamma_trackID"][i])] = int(i)
    al_ev_arr = np.fromiter(ev_g.keys(), dtype=np.int64)
    # accumulate DriftGas daughter hits of Al conv pairs, chunked (bounded memory)
    keep = {}  # pair-index -> {lep: [(time, px,py,pz, arm, ke)]}
    for t in f["HitTree"].iterate(
            ["eventID", "parentID", "particle", "detType", "armID",
             "px", "py", "pz", "time", "ke"], library="np", step_size="250 MB"):
        evc = t["eventID"]
        inev = np.isin(evc, al_ev_arr)          # fast integer prefilter
        if not inev.any():
            continue
        # now only the tiny Al-event subset — string compares are cheap here
        idx = np.where(inev)[0]
        det = t["detType"][idx]
        gasi = (det == b"DriftGas") | (det == "DriftGas")
        if not gasi.any():
            continue
        idx = idx[gasi]
        ev_s = evc[idx]; pid_s = t["parentID"][idx]; part_s = t["particle"][idx]
        for kk, hi in enumerate(idx):
            gmap = ev_g.get(int(ev_s[kk]))
            if not gmap:
                continue
            pi = gmap.get(int(pid_s[kk]))
            if pi is None:
                continue
            p = part_s[kk]
            lep = "e-" if (p == b"e-" or p == "e-") else ("e+" if (p == b"e+" or p == "e+") else None)
            if lep is None:
                continue
            keep.setdefault(pi, {}).setdefault(lep, []).append(
                (t["time"][hi], t["px"][hi], t["py"][hi],
                 t["pz"][hi], t["armID"][hi], t["ke"][hi]))
    for i in al:
        i = int(i)
        rec = {"theta_truth": c["openingAngle"][i],
               "emin_truth": min(c["em_ke"][i], c["ep_ke"][i]),
               "etot_truth": c["em_ke"][i] + c["ep_ke"][i] + 1.022,
               "theta_reco": np.nan, "n_mm": 0, "same_arm": -1,
               "em_arm": -1, "ep_arm": -1, "em_ke_entry": np.nan, "ep_ke_entry": np.nan}
        dirs = {}
        for lep in ("e-", "e+"):
            hits = keep.get(i, {}).get(lep)
            if not hits:
                continue
            arr = np.array(hits)  # cols: time,px,py,pz,arm,ke
            j = int(np.argmin(arr[:, 0]))
            v = arr[j, 1:4]; n = np.linalg.norm(v)
            if n == 0:
                continue
            dirs[lep] = v / n
            if lep == "e-": rec["em_arm"] = int(arr[j, 4]); rec["em_ke_entry"] = float(arr[j, 5])
            else:           rec["ep_arm"] = int(arr[j, 4]); rec["ep_ke_entry"] = float(arr[j, 5])
        rec["n_mm"] = len(dirs)
        if len(dirs) == 2:
            coso = float(np.clip(np.dot(dirs["e-"], dirs["e+"]), -1, 1))
            rec["theta_reco"] = np.degrees(np.arccos(coso))
            rec["same_arm"] = int(rec["em_arm"] == rec["ep_arm"])
        rows.append(rec)
    return len(al)

def main():
    if len(sys.argv) > 1 and sys.argv[1].endswith(".root"):
        files = [sys.argv[1]]; n_events = None
    else:
        files = [f for f in sorted(glob.glob(f"{DIR}/*.root"))]
        # keep only completed files (ConvPairTree present, EventTree full)
        good = []
        n_events = 0
        for f in files:
            try:
                u = uproot.open(f)
                if "ConvPairTree" in [k.split(";")[0] for k in u.keys()] \
                        and u["EventTree"].num_entries >= 9_000_000:
                    good.append(f); n_events += u["EventTree"].num_entries
            except Exception:
                pass
        files = good
        print(f"complete files: {len(files)}  events={n_events:,}")
    rows = []; nal = 0
    for k, fp in enumerate(files):
        nal += process(fp, rows)
        print(f"  [{k+1}/{len(files)}] {fp.split('/')[-1]}  cum Al pairs={len(rows)}", flush=True)
    if not rows:
        print("no Al conv pairs found (campaign not ready?)"); return
    theta_t = np.array([r["theta_truth"] for r in rows])
    theta_r = np.array([r["theta_reco"] for r in rows])
    n_mm    = np.array([r["n_mm"] for r in rows])
    same    = np.array([r["same_arm"] for r in rows])
    emin_t  = np.array([r["emin_truth"] for r in rows])
    pair = n_mm == 2
    print(f"\nAl conv pairs: {len(rows)}")
    print(f"  reached MM (>=1 lepton): {(n_mm>=1).mean()*100:.1f}%")
    print(f"  MM 2-track pair (both):  {pair.mean()*100:.1f}%  (n={pair.sum()})")
    if pair.any():
        print(f"    same-arm (hard to split): {(same[pair]==1).mean()*100:.1f}%")
        print(f"    truth θ (pairs): median={np.median(theta_t[pair]):.1f}  "
              f">108°(X17)={ (theta_t[pair]>108).mean()*100:.2f}%")
        print(f"    RECO  θ (pairs): median={np.median(theta_r[pair]):.1f}  "
              f">108°(X17)={ (theta_r[pair]>108).mean()*100:.2f}%")
        print(f"    reco-truth θ smear: median |Δ|={np.median(np.abs(theta_r[pair]-theta_t[pair])):.1f}°")
    if n_events:
        w = FLUX_WINDOW / n_events
        print(f"\n  RATE (weight/pulse={w:.3e}):")
        print(f"    Al conv pairs produced:  {len(rows)*w:.3g}/pulse  {len(rows)*w*PULSES_DAY:.3g}/day")
        print(f"    Al MM 2-track pairs:     {pair.sum()*w:.3g}/pulse  {pair.sum()*w*PULSES_DAY:.3g}/day")
        if pair.any():
            tail = (pair & (theta_r > 108)).sum()
            print(f"    Al MM pairs reco θ>108° (under X17): {tail*w*PULSES_DAY:.3g}/day")
    np.savez(OUT, theta_t=theta_t, theta_r=theta_r, n_mm=n_mm, same_arm=same,
             emin_t=emin_t,
             em_ke_entry=np.array([r["em_ke_entry"] for r in rows]),
             ep_ke_entry=np.array([r["ep_ke_entry"] for r in rows]),
             n_events=(n_events or 0), flux_window=FLUX_WINDOW, pulses_day=PULSES_DAY)
    print(f"\nwrote {OUT}")

if __name__ == "__main__":
    main()
