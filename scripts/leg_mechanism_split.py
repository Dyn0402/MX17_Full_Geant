#!/usr/bin/env python3
"""leg_mechanism_split.py — classify Al-capture legs (SiPM^plastic per arm,
in-gate, 0.5 MIP) as pair-production vs double-Compton. Signature: a positron
(e+) deposit in the leg arm's SiPM/plastic can only come from pair production;
a leg with electrons only is double Compton. Reports the split and the e+
energy-fraction distribution."""
import glob, sys
import numpy as np, uproot
sys.path.insert(0, "scripts")
from analyze_trigger_thermal import det_class_codes

DIR = sys.argv[1] if len(sys.argv) > 1 else \
      "/eos/experiment/ntof/data/x17/full_sim/neutrons_thermal_trig_2cm"
NF = int(sys.argv[2]) if len(sys.argv) > 2 else 8
FILES = sorted(glob.glob(f"{DIR}/*.root"))[:NF]
MIP_S, MIP_P = 0.458, 4.334 * 2.0 / 2.5
FL, GATE, A, B = 1.41, 1.0, 0.5, 0.5

n_leg = n_pair = 0
epfrac = []
for fp in FILES:
    with uproot.open(fp) as f:
        et = f["EventTree"].arrays(["eventID", "neutron_E_eV", "capture_vol"], library="np")
        cv = np.asarray([x.decode() if isinstance(x, bytes) else x for x in et["capture_vol"]],
                        dtype=object).astype(str)
        mid = int(et["eventID"].max())
        eN = np.zeros(mid + 1); eN[et["eventID"]] = et["neutron_E_eV"]
        isAl = np.zeros(mid + 1, bool); isAl[et["eventID"]] = np.char.startswith(cv, "He3Cap_Al")
        EV = []; ARM = []; CHAN = []; ED = []; EP = []
        for t in f["HitTree"].iterate(["eventID", "armID", "detType", "edep", "particle", "u", "time"],
                                       library="np", step_size="300 MB"):
            code = det_class_codes(t["detType"])
            det = np.where(code == 0, 0, np.where((code == 1) | (code == 2), 1, -1))
            keep = (det >= 0) & isAl[t["eventID"]]
            e_hit = eN[t["eventID"]]
            with np.errstate(divide="ignore"):
                tof = np.where(e_hit > 0, FL / np.sqrt(e_hit), np.inf)
            keep &= (tof + t["time"] * 1e-6) > GATE
            if not keep.any():
                continue
            ev = t["eventID"][keep]; arm = t["armID"][keep]; u = t["u"][keep]
            ed = t["edep"][keep] * 1e-6
            prt = np.asarray([x.decode() if isinstance(x, bytes) else x
                              for x in t["particle"][keep]], dtype=object)
            chan = np.where(det[keep] == 0,
                            np.clip((u + 250.) // 25., 0, 19).astype(int), 20)
            EV.append(ev); ARM.append(arm); CHAN.append(chan); ED.append(ed)
            EP.append((prt == "e+").astype(float) * ed)
        if not EV:
            continue
        ev = np.concatenate(EV); arm = np.concatenate(ARM).astype(np.int64)
        chan = np.concatenate(CHAN).astype(np.int64); ed = np.concatenate(ED); ep = np.concatenate(EP)
        # per (ev,arm,chan) edep -> bar maxima
        bkey = (ev.astype(np.int64) * 4 + arm) * 21 + chan
        ub, bi = np.unique(bkey, return_inverse=True)
        bar = np.zeros(len(ub)); np.add.at(bar, bi, ed)
        akey_b = ub // 21; bchan = ub % 21
        # per (ev,arm) SiPM-max, plastic-max, total edep, e+ edep
        ua, ai_b = np.unique(akey_b, return_inverse=True)
        sipm = np.zeros(len(ua)); plas = np.zeros(len(ua))
        np.maximum.at(sipm, ai_b, np.where(bchan <= 19, bar, 0))
        np.maximum.at(plas, ai_b, np.where(bchan == 20, bar, 0))
        akey_h = ev.astype(np.int64) * 4 + arm
        idx = {k: i for i, k in enumerate(ua)}
        tot = np.zeros(len(ua)); epl = np.zeros(len(ua))
        ai_h = np.array([idx[k] for k in akey_h])
        np.add.at(tot, ai_h, ed); np.add.at(epl, ai_h, ep)
        leg = (sipm / MIP_S >= A) & (plas / MIP_P >= B)
        n_leg += int(leg.sum())
        f_ep = np.where(tot > 0, epl / tot, 0.0)[leg]
        n_pair += int((f_ep > 0.05).sum())
        epfrac.extend(f_ep.tolist())

epfrac = np.array(epfrac)
print(f"files={len(FILES)}  Al legs found: {n_leg}")
print(f"  pair-production legs (e+ energy >5%): {n_pair}  ({n_pair/n_leg*100:.1f}%)")
print(f"  double-Compton legs (e+ <=5%):        {n_leg-n_pair}  ({(n_leg-n_pair)/n_leg*100:.1f}%)")
print(f"  legs with ANY e+ energy:              {(epfrac>0).sum()}  ({(epfrac>0).mean()*100:.1f}%)")
print(f"  median e+ energy fraction (pair legs): {np.median(epfrac[epfrac>0.05])*100:.0f}%")
np.savez("analysis/thermal_2cm/leg_mechanism_split.npz", epfrac=epfrac,
         n_leg=n_leg, n_pair=n_pair)
