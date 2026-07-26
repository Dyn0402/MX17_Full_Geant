#!/usr/bin/env python3
"""
analyze_trigger_provenance.py — what fires the thermal-gate trigger, as a
function of plastic threshold, and how it correlates with the Micromegas.
================================================================================
Trigger (single-arm "leg"): per arm, max single SiPM-bar edep >= 0.5 MIP AND
max plastic-bar edep >= b MIP; the event triggers if >=1 arm fires. b (the
PLASTIC threshold) is scanned; the SiPM threshold is fixed at the nominal
0.5 MIP.

For every triggering event we attribute it to its LEADING firing arm (the arm
with the largest plastic edep among SiPM-passing arms) and trace that arm's
leg-making track with the HitTree birth-truth branches:
  - capture site   : EventTree capture_vol       (what the neutron did)
  - mechanism      : origin_proc                 (compt / conv / ioni / primary)
  - birth volume   : origin_vol   (grouped)      (where the charged track was born)
  - species        : particle
  - MM correlation : does that track cross DriftGas (same-arm), or any DriftGas
                     hit in the event (any-arm)  -> "what we see in the MM"

Sources:
  --thermal  neutrons_thermal_trig_2cm_nose   (gate on; all in-gate; per-pulse rate)
  --epi      neutrons_epi_trig_2cm_nose        (gate on; delayed spill-in only)
  --signal   pairs_thermal_trig_2cm_nose       (no gate; X17=0/IPC=1; efficiency)

Outputs: <outdir>/trigger_provenance.json + trigger_provenance.pdf

Usage (lxplus, LCG_106):
  python3 analyze_trigger_provenance.py \
     --thermal /eos/.../neutrons_thermal_trig_2cm_nose \
     --epi     /eos/.../neutrons_epi_trig_2cm_nose \
     --signal  /eos/.../pairs_thermal_trig_2cm_nose \
     --nthermal 40 --nepi 50 --nsignal 100 --workers 8 -o docs/trigger_provenance
"""
import argparse, glob, json, os, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))
from analyze_trigger_thermal import digest_file  # noqa: E402

# ── calibration + grids (geometry unchanged by the capsule flip) ───────────────
MIP_S, MIP_P = 0.4751, 3.3494   # measured 2 cm muon MPV (2026-07-26 calibration)
SIPM_THR = 0.5                                # fixed SiPM threshold [MIP]
B_GRID = np.round(np.arange(0.1, 3.01, 0.1), 2)   # plastic threshold scan [MIP]
N_PULSE_THERMAL = 4.284e6
N_PULSE_EPI     = 6.484e6

SIPM = ("PlasticScint",)
PLAS = ("BackScintL", "BackScintR")

# capture-site category from EventTree capture_vol
def capture_cat(cv):
    if cv in ("He3Cap_Al", "He3Cap_CFRP"):          return "Al capsule γ"
    if cv.startswith("LiqScint") or cv.startswith("LS_"): return "LS H-capture"
    if cv.startswith("BackScint"):                  return "plastic H-capture"
    if cv == "He3Gas":                              return "gas (n,p)t"
    if cv.startswith(("PCB_", "Micromesh", "Drift", "ResistivePaste",
                      "GasWindow", "AmpGas")):       return "MM / PCB"
    if cv in ("", "World"):                         return "other / none"
    return "other structure"

# birth-volume group from origin_vol (same grouping as the mechanism analysis)
def vol_group(v):
    if v in ("He3Cap_Al", "He3Gas"):    return "Al capsule"
    if v == "He3Cap_CFRP":              return "capsule CFRP"
    if v == "World":                    return "air"
    if v.startswith(("GasWindow", "Drift", "Micromesh", "AmpGas",
                     "ResistivePaste", "PCB_")):   return "MM+PCB"
    if v == "PlasticScint":             return "SiPM bar"
    if v.startswith("BackScint"):       return "plastic+wrap"
    if v.startswith(("LiqScint", "LS_")): return "LS"
    return "other"

CAP_CATS = ["Al capsule γ", "LS H-capture", "plastic H-capture",
            "gas (n,p)t", "MM / PCB", "other structure", "other / none"]


def process_file(fp, gate, want_truth):
    """One file -> per-triggering-candidate-event records (leading firing arm).
    A candidate event has >=1 arm with SiPM >= 0.5 MIP; we record its plastic
    edep [MIP], provenance of the leading arm's leg track, and MM flags."""
    # total generated per event_type (for signal absolute-efficiency denominators)
    ntype = {}
    if want_truth:
        with uproot.open(fp) as f:
            ety = f["EventTree"]["event_type"].array(library="np")
        v, c = np.unique(ety, return_counts=True)
        ntype = {int(k): int(n) for k, n in zip(v, c)}
    obs, aux, nev = digest_file(fp, gate=gate, truth=want_truth)
    if not len(aux["eventID"]):
        return _empty(nev, ntype)
    S = obs[:, :, 0] / MIP_S
    P = obs[:, :, 1] / MIP_P
    sipm_pass = S >= SIPM_THR                      # (N,4)
    cand = sipm_pass.any(axis=1)
    if not cand.any():
        return _empty(nev, ntype)

    # leading firing arm = argmax plastic edep among SiPM-passing arms
    Pmask = np.where(sipm_pass, P, -1.0)
    arm_lead = Pmask.argmax(axis=1)
    P_trig = Pmask[np.arange(len(P)), arm_lead]
    ev_c   = aux["eventID"][cand]
    arm_c  = arm_lead[cand]
    Ptrig_c = P_trig[cand]
    etype_c = (aux["event_type"][cand] if want_truth else
               np.full(cand.sum(), 2, np.int32))

    # re-read HitTree only for candidate events; classify the leading arm
    ev_set = set(int(e) for e in ev_c)
    with uproot.open(fp) as f:
        et = f["EventTree"].arrays(["eventID", "capture_vol"], library="np")
        capmap = {int(e): str(v).rstrip("\x00")
                  for e, v in zip(et["eventID"], et["capture_vol"])
                  if int(e) in ev_set}
        cols = ["eventID", "armID", "detType", "trackID", "parentID",
                "particle", "edep", "origin_vol", "origin_proc", "origin_ke", "ke"]
        REL_DET = list(SIPM) + list(PLAS) + ["DriftGas"]
        parts = []
        for t in f["HitTree"].iterate(cols, library="np", step_size="80 MB"):
            m = np.isin(t["eventID"], ev_c)          # cheap int filter first
            if not m.any():
                continue
            sub = {k: v[m] for k, v in t.items()}
            det_c = np.array([str(x).rstrip("\x00") for x in sub["detType"]])
            code = np.full(len(det_c), -1, np.int8)  # 0=SiPM 1=plastic 2=DriftGas
            code[det_c == "PlasticScint"] = 0
            code[np.isin(det_c, ("BackScintL", "BackScintR"))] = 1
            code[det_c == "DriftGas"] = 2
            rel = code >= 0                           # keep only SiPM/plastic/MM hits
            if rel.any():
                del sub["detType"]                   # drop the big string column
                d = {k: v[rel] for k, v in sub.items()}
                d["detcode"] = code[rel]
                parts.append(d)
    if not parts:
        return _empty(nev, ntype)
    h = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    hev, harm = h["eventID"], h["armID"]
    dc = h["detcode"]
    is_dg = dc == 2
    is_s = dc == 0
    is_p = dc == 1

    rec = {k: [] for k in ("P_trig", "etype", "capcat", "proc", "birthvol",
                           "species", "birth_ke", "mm_same", "mm_any", "same_track")}
    # group hit indices by event for speed
    order = np.argsort(hev, kind="stable")
    hev_s = hev[order]
    lo_all = np.searchsorted(hev_s, ev_c, "left")
    hi_all = np.searchsorted(hev_s, ev_c, "right")

    for k, (ev, arm, pt, ety) in enumerate(zip(ev_c, arm_c, Ptrig_c, etype_c)):
        idx = order[lo_all[k]:hi_all[k]]           # all hits for this event
        a_idx = idx[harm[idx] == arm]              # leading-arm hits
        ss = a_idx[is_s[a_idx]]
        pp = a_idx[is_p[a_idx]]
        ts = set(h["trackID"][ss].tolist()); tp = set(h["trackID"][pp].tolist())
        common = ts & tp
        # leg track: shared SiPM&plastic track (>=98%), else the top-edep plastic track
        if common:
            same = 1
            # pick common track with largest plastic deposit
            best, bde = None, -1.0
            for tid in common:
                de = h["edep"][pp[h["trackID"][pp] == tid]].sum()
                if de > bde: bde, best = de, tid
            legtrk = best
            j = ss[h["trackID"][ss] == legtrk][0]  # SiPM hit of the leg track
        else:
            same = 0
            if len(pp) == 0:
                continue
            legtrk = int(h["trackID"][pp[np.argmax(h["edep"][pp])]])
            j = pp[h["trackID"][pp] == legtrk][0]
        rec["P_trig"].append(float(pt))
        rec["etype"].append(int(ety))
        rec["capcat"].append(capture_cat(capmap.get(int(ev), "")))
        rec["proc"].append(str(h["origin_proc"][j]).rstrip("\x00"))
        rec["birthvol"].append(vol_group(str(h["origin_vol"][j]).rstrip("\x00")))
        rec["species"].append(str(h["particle"][j]).rstrip("\x00"))
        rec["birth_ke"].append(float(h["origin_ke"][j]))
        rec["same_track"].append(same)
        # MM: leg track crosses drift gas (same-arm) / any drift-gas hit (any-arm)
        dg_ev = idx[is_dg[idx]]
        rec["mm_any"].append(int(len(dg_ev) > 0))
        rec["mm_same"].append(int((h["trackID"][dg_ev] == legtrk).any()
                                  if len(dg_ev) else False))
    out = {k: np.array(v) for k, v in rec.items()}
    out["nev"] = nev
    out["ntype"] = ntype
    return out


def _empty(nev, ntype=None):
    d = {k: np.array([]) for k in ("P_trig", "etype", "capcat", "proc",
         "birthvol", "species", "birth_ke", "mm_same", "mm_any", "same_track")}
    d["nev"] = nev
    d["ntype"] = ntype or {}
    return d


def run_source(files, gate, want_truth, workers, label):
    agg = {k: [] for k in ("P_trig", "etype", "capcat", "proc", "birthvol",
                           "species", "birth_ke", "mm_same", "mm_any", "same_track")}
    nev_tot = 0
    ntype_tot = {}
    def absorb(r):
        nonlocal nev_tot
        nev_tot += r["nev"]
        for k in agg: agg[k].append(r[k])
        for t, n in r.get("ntype", {}).items():
            ntype_tot[t] = ntype_tot.get(t, 0) + n
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(process_file, fp, gate, want_truth): fp for fp in files}
            for i, fut in enumerate(as_completed(futs)):
                absorb(fut.result())
                print(f"  [{label}] {i+1}/{len(files)} cum_cand="
                      f"{sum(len(a) for a in agg['P_trig'])}", flush=True)
    else:
        for i, fp in enumerate(files):
            absorb(process_file(fp, gate, want_truth))
            print(f"  [{label}] {i+1}/{len(files)} cum_cand="
                  f"{sum(len(a) for a in agg['P_trig'])}", flush=True)
    out = {k: (np.concatenate(v) if v else np.array([])) for k, v in agg.items()}
    out["nev_tot"] = nev_tot
    out["ntype"] = ntype_tot
    return out


def collect(d, n):
    fs = sorted(glob.glob(f"{d}/*.root"))
    return fs[:n] if n else fs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--thermal"); ap.add_argument("--epi"); ap.add_argument("--signal")
    ap.add_argument("--nthermal", type=int, default=0)
    ap.add_argument("--nepi", type=int, default=0)
    ap.add_argument("--nsignal", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("-o", "--outdir", default="trigger_provenance")
    args = ap.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

    res = {}
    if args.thermal:
        print("THERMAL:"); res["thermal"] = run_source(
            collect(args.thermal, args.nthermal), True, False, args.workers, "th")
    if args.epi:
        print("EPI:"); res["epi"] = run_source(
            collect(args.epi, args.nepi), True, False, args.workers, "epi")
    if args.signal:
        print("SIGNAL:"); res["signal"] = run_source(
            collect(args.signal, args.nsignal), False, True, args.workers, "sig")

    summary = build_summary(res)
    json.dump(summary, open(out / "trigger_provenance.json", "w"), indent=1,
              default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o)
    print("wrote", out / "trigger_provenance.json")
    try:
        make_figures(res, summary, out / "trigger_provenance.pdf")
        print("wrote", out / "trigger_provenance.pdf")
    except Exception as e:
        print("figure step failed:", e)


def rate_vs_b(P_trig, npulse, nev_tot, mask=None):
    p = P_trig if mask is None else P_trig[mask]
    return np.array([(p >= b).sum() for b in B_GRID]) / nev_tot * npulse


def build_summary(res):
    s = {"B_GRID": B_GRID.tolist(), "MIP_S": MIP_S, "MIP_P": MIP_P,
         "SIPM_THR": SIPM_THR, "b_nom": 0.5}
    inom = int(np.argmin(np.abs(B_GRID - 0.5)))
    if "thermal" in res:
        th = res["thermal"]; nev = th["nev_tot"]
        s["thermal"] = {"nev_tot": nev, "n_cand": int(len(th["P_trig"]))}
        s["thermal"]["rate_total"] = rate_vs_b(th["P_trig"], N_PULSE_THERMAL, nev).tolist()
        # per capture-site category
        s["thermal"]["rate_by_capcat"] = {
            c: rate_vs_b(th["P_trig"], N_PULSE_THERMAL, nev, th["capcat"] == c).tolist()
            for c in CAP_CATS}
        # MM-tagged rates vs b (same-arm / any-arm)
        s["thermal"]["rate_mm_same"] = rate_vs_b(th["P_trig"], N_PULSE_THERMAL, nev,
                                                 th["mm_same"] == 1).tolist()
        s["thermal"]["rate_mm_any"] = rate_vs_b(th["P_trig"], N_PULSE_THERMAL, nev,
                                                th["mm_any"] == 1).tolist()
        # provenance at nominal threshold b=0.5
        trig = th["P_trig"] >= 0.5
        s["thermal"]["rate_nom"] = float(trig.sum() / nev * N_PULSE_THERMAL)
        def tally(field):
            vals, cnts = np.unique(th[field][trig], return_counts=True)
            order = np.argsort(-cnts)
            return {str(vals[i]): int(cnts[i]) for i in order}
        s["thermal"]["nom_capcat"] = tally("capcat")
        s["thermal"]["nom_proc"] = tally("proc")
        s["thermal"]["nom_birthvol"] = tally("birthvol")
        s["thermal"]["nom_species"] = tally("species")
        s["thermal"]["nom_ntrig"] = int(trig.sum())
        s["thermal"]["nom_mm_same_frac"] = float(th["mm_same"][trig].mean())
        s["thermal"]["nom_mm_any_frac"] = float(th["mm_any"][trig].mean())
        s["thermal"]["nom_same_track_frac"] = float(th["same_track"][trig].mean())
        # MM-tag fraction vs b (any-arm) for the total
        with np.errstate(invalid="ignore"):
            s["thermal"]["mm_any_frac_vs_b"] = [
                float(th["mm_any"][th["P_trig"] >= b].mean())
                if (th["P_trig"] >= b).any() else float("nan") for b in B_GRID]
            s["thermal"]["mm_same_frac_vs_b"] = [
                float(th["mm_same"][th["P_trig"] >= b].mean())
                if (th["P_trig"] >= b).any() else float("nan") for b in B_GRID]
    if "epi" in res:
        ep = res["epi"]; nev = ep["nev_tot"]
        s["epi"] = {"nev_tot": nev, "n_cand": int(len(ep["P_trig"])),
                    "rate_total": rate_vs_b(ep["P_trig"], N_PULSE_EPI, nev).tolist(),
                    "rate_nom": float((ep["P_trig"] >= 0.5).sum() / nev * N_PULSE_EPI)}
    if "signal" in res:
        sg = res["signal"]; ntype = sg.get("ntype", {})
        for lbl, code in (("X17", 0), ("IPC", 1)):
            m = sg["etype"] == code
            ng = int(ntype.get(code, ntype.get(str(code), m.sum())))  # total GENERATED
            eff = np.array([(sg["P_trig"][m] >= b).sum() for b in B_GRID]) / max(ng, 1)
            trig = m & (sg["P_trig"] >= 0.5)
            s.setdefault("signal", {})[lbl] = {
                "n_gen": ng, "n_cand": int(m.sum()),
                "eff_vs_b": eff.tolist(),
                "mm_any_frac_nom": float(sg["mm_any"][trig].mean()) if trig.any() else None,
                "mm_same_frac_nom": float(sg["mm_same"][trig].mean()) if trig.any() else None,
            }
    return s


def make_figures(res, s, pdf_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    b = np.array(s["B_GRID"]); bMeV = b * MIP_P
    COL = {"Al capsule γ": "#e34948", "LS H-capture": "#2a78d6",
           "plastic H-capture": "#1baf7a", "gas (n,p)t": "#eb6834",
           "MM / PCB": "#8a6d3b", "other structure": "#9b59b6",
           "other / none": "#6b6a66"}
    with PdfPages(pdf_path) as pdf:
        # PAGE 1: stacked rate vs plastic threshold + signal efficiency
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
        th = s.get("thermal", {})
        if th:
            cats = [c for c in CAP_CATS if any(np.array(th["rate_by_capcat"][c]) > 0)]
            ys = np.vstack([th["rate_by_capcat"][c] for c in cats])
            ax1.stackplot(b, ys, labels=cats, colors=[COL[c] for c in cats], alpha=0.9)
            ax1.plot(b, th["rate_total"], "k-", lw=1.4, label="total")
            if "epi" in s:
                ax1.plot(b, s["epi"]["rate_total"], "k--", lw=1.1, label="epi spill-in")
        ax1.set_yscale("log"); ax1.set_ylim(0.05, max(th.get("rate_total", [1]))*1.5)
        ax1.set_xlabel("plastic threshold b  [MIP]")
        ax1.set_ylabel("trigger rate  [legs / pulse, in-gate]")
        ax1.set_title("What fires the trigger vs plastic threshold\n(single-arm SiPM∧plastic, SiPM≥0.5 MIP)")
        ax1.legend(fontsize=7.5, loc="upper right")
        secax = ax1.secondary_xaxis("top", functions=(lambda x: x*MIP_P, lambda x: x/MIP_P))
        secax.set_xlabel("plastic threshold  [MeV]")
        ax1.axvline(0.5, color="0.5", ls=":", lw=0.8)
        # signal efficiency
        if "signal" in s:
            for lbl, c in (("X17", "#e34948"), ("IPC", "#2a78d6")):
                ax2.plot(b, s["signal"][lbl]["eff_vs_b"], color=c, lw=1.6, label=f"{lbl} ε")
        ax2.set_xlabel("plastic threshold b  [MIP]"); ax2.set_ylabel("signal trigger efficiency")
        ax2.set_title("Signal survival vs plastic threshold"); ax2.set_ylim(0, 1)
        ax2.axvline(0.5, color="0.5", ls=":", lw=0.8); ax2.legend(fontsize=9)
        secax2 = ax2.secondary_xaxis("top", functions=(lambda x: x*MIP_P, lambda x: x/MIP_P))
        secax2.set_xlabel("plastic threshold  [MeV]")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # PAGE 2: provenance breakdown at nominal threshold
        if th:
            fig, axs = plt.subplots(2, 2, figsize=(13, 8))
            for ax, key, title in ((axs[0,0], "nom_capcat", "capture site (what the neutron did)"),
                                   (axs[0,1], "nom_proc", "birth process (mechanism)"),
                                   (axs[1,0], "nom_birthvol", "birth volume of the leg track"),
                                   (axs[1,1], "nom_species", "species crossing the trigger")):
                d = th[key]; tot = sum(d.values())
                items = list(d.items())[:9]
                labs = [k for k, _ in items]; vals = [v/tot*100 for _, v in items]
                y = np.arange(len(labs))[::-1]
                ax.barh(y, vals, color="#2a78d6", height=0.6)
                for yy, v in zip(y, vals): ax.text(v+0.5, yy, f"{v:.1f}%", va="center", fontsize=8)
                ax.set_yticks(y, labs, fontsize=8.5); ax.set_xlim(0, max(vals)*1.18)
                ax.set_title(title, fontsize=10); ax.spines[["top","right"]].set_visible(False)
            fig.suptitle(f"How the trigger background gets made — nominal b=0.5 MIP "
                         f"({th['rate_nom']:.0f} legs/pulse; {th['nom_same_track_frac']*100:.0f}% one crossing track)",
                         fontsize=12)
            fig.tight_layout(rect=(0,0,1,0.96)); pdf.savefig(fig); plt.close(fig)

        # PAGE 3: MM correlation — trigger vs Micromegas
        if th:
            fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
            axL.plot(b, th["rate_total"], "k-", lw=1.5, label="all triggers")
            axL.plot(b, th["rate_mm_any"], color="#1baf7a", lw=1.5, label="+ MM track (any-arm)")
            axL.plot(b, th["rate_mm_same"], color="#e34948", lw=1.5, label="+ MM track (trigger arm)")
            axL.fill_between(b, th["rate_mm_same"], th["rate_total"], color="0.85",
                             label="MM-blind (no track on trigger arm)")
            axL.set_yscale("log"); axL.set_xlabel("plastic threshold b  [MIP]")
            axL.set_ylabel("rate [legs / pulse]"); axL.legend(fontsize=8)
            axL.set_title("Trigger vs Micromegas: how much is MM-tagged")
            axL.axvline(0.5, color="0.5", ls=":", lw=0.8)
            frac = np.array(th["mm_any_frac_vs_b"]); fracs = np.array(th["mm_same_frac_vs_b"])
            axR.plot(b, frac*100, color="#1baf7a", lw=1.6, label="MM any-arm")
            axR.plot(b, fracs*100, color="#e34948", lw=1.6, label="MM trigger-arm")
            if "signal" in s:
                for lbl, c in (("X17", "#111111"),):
                    v = s["signal"][lbl].get("mm_any_frac_nom")
                    if v is not None:
                        axR.axhline(v*100, color=c, ls="--", lw=1.0,
                                    label=f"{lbl} signal (MM any, b=0.5)")
            axR.set_xlabel("plastic threshold b  [MIP]")
            axR.set_ylabel("% of triggers with an MM track"); axR.set_ylim(0, 100)
            axR.set_title("MM-tag fraction of triggers vs threshold"); axR.legend(fontsize=8)
            axR.axvline(0.5, color="0.5", ls=":", lw=0.8)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


if __name__ == "__main__":
    main()
