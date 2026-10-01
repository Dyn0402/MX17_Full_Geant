#!/usr/bin/env python3
"""
ill_pairs.py -- opening-angle reconstruction and acceptance for the ILL scan
============================================================================
HANDOFF_SIM.md §5g/§7 (x17_facility_search/ill).  Two stages:

    # per ROOT file (condor) -> compact per-event table (.npz)
    python3 ill_pairs.py reduce <pairs_jobNNN_t0.root> -o <part.npz>

    # merge -> per-estimator σ68 / bias and acceptance × ε, JSON + CSV
    python3 ill_pairs.py merge <parts/*.npz> -o <outdir> --config G1 \\
        [--assumed-vertex 0,0,0] [--dir-smear-deg 10]

reduce keeps, per event: truth (type, mass, opening angle, lepton KE and
direction, vertex) and, per lepton, its FIRST DriftGas hit (position, and the
momentum direction there), the edep-weighted centroid and a PCA line fit of
its drift hits in that arm (positions smeared by 0.5 mm, as analyze_pairs.py),
plus the arm bitmask of trigger legs (SiPM bar ≥ 0.5 MIP AND plastic bar
≥ 0.5 MIP in one arm, thermal_accounting definitions).

merge evaluates every estimator on the selected events (both leptons in the
gaps, in different arms, and a pair-tag: legs in ≥ 2 arms):

    vline      chord to the TRUE vertex (oracle)
    nomline    chord to an assumed vertex (default (0,0,0): beam axis, y = 0)
    first      momentum direction at the first gap hit (direction only)
    fit        PCA line of the smeared drift hits (direction only)
    xing       chord to (0, y*, 0), y* = mean of the two fitted lines'
               closest-approach y to the beam axis (per-event vertex)
    xing_s     the same with each fitted direction smeared by --dir-smear-deg
               (a data-like direction error; the simulation's PCA direction is
               far better than the chambers' measured 11–26°)

σ68 = half the 16–84 % width of (θ_reco − θ_truth); bias = its median.
Reported overall, in the 100–180° truth window and against the softer
lepton's kinetic energy.
"""
from __future__ import annotations

import argparse
import glob
import json
import zlib
from pathlib import Path

import numpy as np

MM_HIT_SMEAR_MM = 0.5
MIP_S_EV, MIP_P_EV, LEG_MIP = 0.4751e6, 3.3494e6, 0.5
KE_BINS = np.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 21.0])
EST = ("vline", "nomline", "first", "fit", "xing", "xing_s")
BRANCHES = ["eventID", "trackID", "parentID", "armID", "detType", "particle",
            "u", "edep", "gx", "gy", "gz", "px", "py", "pz"]


def _s(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


# --------------------------------------------------------------------------- #
# reduce
# --------------------------------------------------------------------------- #
def _lepton_table(ev, arm, X, D, edep, rng):
    """Per event (one lepton): first-hit arm/position/direction, and over that
    arm's hits the edep centroid and the PCA fit of 0.5 mm-smeared positions.
    Inputs are that lepton's DriftGas hits in stepping order."""
    ue, first = np.unique(ev, return_index=True)
    out = dict(ev=ue, arm=arm[first], P=X[first], d0=D[first])
    keep = arm == np.repeat(out["arm"], np.diff(np.r_[first, len(ev)]))
    ev, X, edep, D = ev[keep], X[keep], edep[keep], D[keep]
    _, inv = np.unique(ev, return_inverse=True)
    n = np.bincount(inv, minlength=len(ue)).astype(float)
    w = np.bincount(inv, edep, minlength=len(ue)).clip(1e-30)
    out["C"] = np.stack([np.bincount(inv, edep * X[:, k], minlength=len(ue)) for k in range(3)], 1) / w[:, None]
    Xs = X + rng.normal(0, MM_HIT_SMEAR_MM, X.shape)
    mean = np.stack([np.bincount(inv, Xs[:, k], minlength=len(ue)) for k in range(3)], 1) / n[:, None]
    Xc = Xs - mean[inv]
    M = np.empty((len(ue), 3, 3))
    for a in range(3):
        for b in range(a, 3):
            M[:, a, b] = M[:, b, a] = np.bincount(inv, Xc[:, a] * Xc[:, b], minlength=len(ue))
    _, vecs = np.linalg.eigh(M)
    d = vecs[:, :, 2]
    sgn = np.sign((d * out["d0"]).sum(1))
    sgn[sgn == 0] = 1
    d = d * sgn[:, None]
    d[n < 3] = np.nan
    out["fit"], out["Cs"], out["n"] = d, mean, n
    return out


def reduce_file(fp: str, out: str, step="200 MB"):
    import uproot
    rng = np.random.default_rng(zlib.crc32(Path(fp).name.encode()))
    with uproot.open(fp) as f:
        et = f["EventTree"].arrays(["eventID", "event_type", "inv_mass", "openingAngle",
                                    "em_ke", "ep_ke", "em_px", "em_py", "em_pz",
                                    "ep_px", "ep_py", "ep_pz", "vtx_x", "vtx_y", "vtx_z"],
                                   library="np")
        N = len(et["eventID"])
        emax = int(et["eventID"].max()) + 1
        lep = {p: dict(arm=np.full(emax, -1, np.int8), P=np.full((emax, 3), np.nan, np.float32),
                       d0=np.full((emax, 3), np.nan, np.float32),
                       C=np.full((emax, 3), np.nan, np.float32),
                       Cs=np.full((emax, 3), np.nan, np.float32),
                       fit=np.full((emax, 3), np.nan, np.float32),
                       n=np.zeros(emax, np.int32)) for p in ("em", "ep")}
        sums = {}                     # (event*4+arm)*22+chan -> edep

        def process(t):
            det = _s(t["detType"])
            ev = t["eventID"].astype(np.int64)
            arm = t["armID"].astype(np.int64)
            # scintillator channel sums for legs (all particles)
            sc = (det == "PlasticScint") | (det == "BackScintL") | (det == "BackScintR")
            if sc.any():
                ch = np.where(det[sc] == "PlasticScint",
                              np.clip((t["u"][sc] + 250.0) // 25.0, 0, 19),
                              np.where(det[sc] == "BackScintL", 20, 21)).astype(np.int64)
                k = (ev[sc] * 4 + arm[sc]) * 22 + ch
                uk, inv = np.unique(k, return_inverse=True)
                s = np.bincount(inv, t["edep"][sc])
                for kk, ss in zip(uk.tolist(), s.tolist()):
                    sums[kk] = sums.get(kk, 0.0) + ss
            # primary leptons in the drift gaps
            prt = _s(t["particle"])
            for p, name in (("em", "e-"), ("ep", "e+")):
                m = (det == "DriftGas") & (t["parentID"] == 0) & (prt == name)
                if not m.any():
                    continue
                X = np.stack([t["gx"][m], t["gy"][m], t["gz"][m]], 1).astype(float)
                D = np.stack([t["px"][m], t["py"][m], t["pz"][m]], 1).astype(float)
                L = _lepton_table(ev[m], arm[m], X, D, t["edep"][m].astype(float), rng)
                for key in ("arm", "P", "d0", "C", "Cs", "fit", "n"):
                    lep[p][key][L["ev"]] = L[key]

        # hits are written event by event (single-threaded), so eventID is
        # non-decreasing: hold back the last event of each chunk until the next
        carry = None
        for t in f["HitTree"].iterate(BRANCHES, library="np", step_size=step):
            if carry is not None:
                t = {k: np.concatenate([carry[k], t[k]]) for k in t}
            cut = np.searchsorted(t["eventID"], t["eventID"][-1], side="left")
            carry = {k: v[cut:] for k, v in t.items()}
            if cut:
                process({k: v[:cut] for k, v in t.items()})
        if carry is not None and len(carry["eventID"]):
            process(carry)
    # legs
    legmask = np.zeros(emax, np.int8)
    if sums:
        k = np.fromiter(sums.keys(), np.int64)
        s = np.fromiter(sums.values(), float)
        s_ok = (k % 22 < 20) & (s >= LEG_MIP * MIP_S_EV)
        p_ok = (k % 22 >= 20) & (s >= LEG_MIP * MIP_P_EV)
        ea = np.intersect1d(np.unique(k[s_ok] // 22), np.unique(k[p_ok] // 22))
        np.bitwise_or.at(legmask, ea // 4, (1 << (ea % 4)).astype(np.int8))
    e = et["eventID"].astype(np.int64)
    rec = dict(eventID=e, event_type=et["event_type"].astype(np.int8),
               inv_mass=et["inv_mass"].astype(np.float32),
               theta=et["openingAngle"].astype(np.float32),
               em_ke=et["em_ke"].astype(np.float32), ep_ke=et["ep_ke"].astype(np.float32),
               em_d=np.stack([et["em_px"], et["em_py"], et["em_pz"]], 1).astype(np.float32),
               ep_d=np.stack([et["ep_px"], et["ep_py"], et["ep_pz"]], 1).astype(np.float32),
               V=np.stack([et["vtx_x"], et["vtx_y"], et["vtx_z"]], 1).astype(np.float32),
               legs=legmask[e])
    for p in ("em", "ep"):
        for key, v in lep[p].items():
            rec[f"{p}_{key}"] = v[e]
    np.savez_compressed(out, **rec)
    print(f"{fp}: {N} events -> {out}")


# --------------------------------------------------------------------------- #
# merge
# --------------------------------------------------------------------------- #
def _ang(a, b):
    na = np.linalg.norm(a, axis=1)
    nb = np.linalg.norm(b, axis=1)
    c = (a * b).sum(1) / (na * nb)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def _smear_dir(d, deg, rng):
    """Rotate each unit vector by a Gaussian polar angle (σ = deg), random azimuth."""
    if deg <= 0:
        return d
    d = d / np.linalg.norm(d, axis=1)[:, None]
    a = np.where(np.abs(d[:, 0]) < 0.9, 1.0, 0.0)
    ref = np.stack([a, 1 - a, np.zeros_like(a)], 1)
    e1 = np.cross(d, ref)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(d, e1)
    th = np.radians(rng.normal(0, deg, len(d)))
    ph = rng.uniform(0, 2 * np.pi, len(d))
    return (np.cos(th)[:, None] * d + np.sin(th)[:, None] *
            (np.cos(ph)[:, None] * e1 + np.sin(ph)[:, None] * e2))


def _y_star(C, d):
    """y of each line's closest approach to the beam axis (x = z = 0)."""
    den = d[:, 0] ** 2 + d[:, 2] ** 2
    t = -(C[:, 0] * d[:, 0] + C[:, 2] * d[:, 2]) / np.where(den > 0, den, np.nan)
    return C[:, 1] + t * d[:, 1]


def _stats(delta):
    delta = delta[np.isfinite(delta)]
    if len(delta) < 20:
        return dict(n=int(len(delta)), sigma68=None, bias=None, rms=None)
    q16, q50, q84 = np.percentile(delta, [16, 50, 84])
    return dict(n=int(len(delta)), sigma68=float(0.5 * (q84 - q16)), bias=float(q50),
                rms=float(np.sqrt(np.mean(delta ** 2))))


def merge(parts, outdir, config, assumed, dir_smear, seed=5):
    R = {}
    for p in parts:
        z = np.load(p)
        for k in z.files:
            R.setdefault(k, []).append(z[k])
    R = {k: np.concatenate(v) for k, v in R.items()}
    rng = np.random.default_rng(seed)
    N = len(R["eventID"])
    smear = lambda X: X + rng.normal(0, MM_HIT_SMEAR_MM, X.shape)   # noqa: E731
    Pm, Pp = smear(R["em_P"].astype(float)), smear(R["ep_P"].astype(float))
    V = R["V"].astype(float)
    A = np.array(assumed, float)[None, :]
    fm, fp = R["em_fit"].astype(float), R["ep_fit"].astype(float)
    Cm, Cp = R["em_Cs"].astype(float), R["ep_Cs"].astype(float)
    ym, yp = _y_star(Cm, fm), _y_star(Cp, fp)
    ys = 0.5 * (ym + yp)
    fms, fps = _smear_dir(fm, dir_smear, rng), _smear_dir(fp, dir_smear, rng)
    ys_s = 0.5 * (_y_star(Cm, fms) + _y_star(Cp, fps))
    Vx = np.stack([np.zeros(N), ys, np.zeros(N)], 1)
    Vxs = np.stack([np.zeros(N), ys_s, np.zeros(N)], 1)
    reco = dict(vline=_ang(Pm - V, Pp - V), nomline=_ang(Pm - A, Pp - A),
                first=_ang(R["em_d0"], R["ep_d0"]), fit=_ang(fm, fp),
                xing=_ang(Pm - Vx, Pp - Vx), xing_s=_ang(Pm - Vxs, Pp - Vxs))
    th = R["theta"].astype(float)
    mm2 = (R["em_arm"] >= 0) & (R["ep_arm"] >= 0)
    arm2 = mm2 & (R["em_arm"] != R["ep_arm"])
    nleg = np.array([bin(int(x) & 0xF).count("1") for x in R["legs"]])
    tag = nleg >= 2
    sel = arm2 & tag
    softer = np.minimum(R["em_ke"], R["ep_ke"]).astype(float)
    win = (th >= 100) & (th <= 180)

    out = dict(config=config, assumed_vertex=list(assumed), dir_smear_deg=dir_smear,
               hit_smear_mm=MM_HIT_SMEAR_MM, channels={})
    names = {0: "X17", 1: "IPC"}
    for et in sorted(set(R["event_type"].tolist())):
        m = R["event_type"] == et
        n = int(m.sum())
        ch = dict(n_generated=n,
                  acc_mm_both=float((mm2 & m).sum() / n), acc_mm_2arm=float((arm2 & m).sum() / n),
                  eff_pairtag=float((tag & m).sum() / n), acc_x_eff=float((sel & m).sum() / n),
                  acc_x_eff_100_180=float((sel & m & win).sum() / max(1, (m & win).sum())),
                  frac_truth_100_180=float((m & win).sum() / n),
                  vertex_sigma_along_mm=float(V[m, 1].std()),
                  vertex_sigma_across_mm=float(np.sqrt(0.5 * (V[m, 0].var() + V[m, 2].var()))),
                  estimators={})
        dy = ys[sel & m] - V[sel & m, 1]
        ch["xing_vertex_y"] = _stats(dy)
        for e in EST:
            d = reco[e] - th
            row = dict(all=_stats(d[sel & m]), win_100_180=_stats(d[sel & m & win]), vs_softer_ke=[])
            for lo, hi in zip(KE_BINS[:-1], KE_BINS[1:]):
                mk = sel & m & (softer >= lo) & (softer < hi)
                row["vs_softer_ke"].append(dict(ke_lo=float(lo), ke_hi=float(hi), **_stats(d[mk])))
            ch["estimators"][e] = row
        out["channels"][names.get(et, str(et))] = ch
    od = Path(outdir)
    od.mkdir(parents=True, exist_ok=True)
    (od / "pairs.json").write_text(json.dumps(out, indent=1))
    print(f"{config}: {N:,} events, assumed vertex {assumed}, dir smear {dir_smear}°")
    for chn, ch in out["channels"].items():
        print(f"  {chn}: n={ch['n_generated']:,}  MM both {ch['acc_mm_both']:.3f}  2-arm "
              f"{ch['acc_mm_2arm']:.3f}  pair-tag {ch['eff_pairtag']:.3f}  acc×ε {ch['acc_x_eff']:.4f}"
              f"  (100–180°: {ch['acc_x_eff_100_180']:.4f})  vtx σ along/across "
              f"{ch['vertex_sigma_along_mm']:.1f}/{ch['vertex_sigma_across_mm']:.1f} mm")
        for e in EST:
            a, w = ch["estimators"][e]["all"], ch["estimators"][e]["win_100_180"]
            f = lambda s: "   n/a" if s["sigma68"] is None else f"{s['sigma68']:6.2f}/{s['bias']:+6.2f}"  # noqa: E731
            print(f"     {e:8s} σ68/bias all {f(a)}   100–180° {f(w)}  (n={w['n']})")
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
    m.add_argument("--assumed-vertex", default="0,0,0")
    m.add_argument("--dir-smear-deg", type=float, default=10.0)
    a = ap.parse_args()
    if a.cmd == "reduce":
        reduce_file(a.file, a.out)
    else:
        files = []
        for p in a.parts:
            files += sorted(glob.glob(p)) if "*" in p else [p]
        merge(files, a.outdir, a.config, [float(x) for x in a.assumed_vertex.split(",")],
              a.dir_smear_deg)


if __name__ == "__main__":
    main()
