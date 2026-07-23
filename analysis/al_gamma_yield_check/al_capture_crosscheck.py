#!/usr/bin/env python3
"""Analytic cross-check of the 27Al(n,gamma) capture rate in the He-3 capsule
aluminium pressure vessel vs Geant4 (.claude/al_gamma_yield_check/HANDOFF.md).

Single-pass (straight-line attenuation) transport folded over
  - the exact STEP-derived polycone geometry mirrored from
    src/DetectorConstruction.cc:341-412, in the AS-BUILT orientation: the
    placement rotation rotateX(-90 deg) maps local z -> world -y, so the
    VALVE (local z=+51) faces the incoming beam, not the nose;
  - the EAR2 flux histogram (energy) x Lambda2D (radius | energy), sampled
    exactly as src/X17PrimaryGenerator.cc does (full weight of window-
    overlapping bins, nearest profile row in log10 E, 0.5 mm radial bins);
  - 1/v capture cross section for 27Al anchored at sigma0 = 0.231 b
    (2200 m/s; consistent with the 1778.9 keV decay-line sigma_gamma =
    0.232(4) b of the IAEA PGAA k0 table);
  - He-3 (n,p) removal from ENDF via data/He3.h5 (MT=103).

Outputs docs/al_gamma_yield_check/analytic_results.npz + printed tables.
"""
from pathlib import Path
import json
import numpy as np
import h5py
import uproot

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent                 # repo root (analysis/al_gamma_yield_check/..)
OUT = ROOT / "docs/al_gamma_yield_check"
OUT.mkdir(parents=True, exist_ok=True)

# ── constants ────────────────────────────────────────────────────────────────
NA = 6.02214076e23
RHO_AL, A_AL = 2.699, 26.9815385          # g/cm3, g/mol  (G4_Al)
RHO_HE3, A_HE3 = 0.0627, 3.0160293       # 500 bar fill
N_AL = RHO_AL * NA / A_AL                 # 6.024e22 cm^-3
N_HE3 = RHO_HE3 * NA / A_HE3              # 1.252e22 cm^-3
SIG0_AL = 0.231e-24                       # cm2 @ 0.0253 eV  (Mughabghab / EGAF)
E0 = 0.0253                               # eV
I_7724 = 0.0493 / 0.231                   # 7.724 MeV gammas per capture (EGAF)
PULSES_DAY = 1.929e4

# ── geometry: polycone vertices [mm], local frame (z=-35 tip .. +51 valve) ──
zGas = np.array([-29.5, -28, -26, -24, -22, -20, -15, -5, 5, 15, 20, 22, 24,
                 26, 28, 30, 32, 34, 36, 38, 40, 44, 50.7])
roGas = np.array([0.001, 6.0, 8.0, 9.165, 9.798, 10.0, 10.0, 10.0, 10.0, 10.0,
                  10.0, 9.798, 9.165, 8.0, 6.299, 4.842, 3.660, 2.711, 1.967,
                  1.410, 1.026, 0.750, 0.750])
zAl = np.array([-35, -34, -33, -31, -29, -27, -25, -23, -21, -20, -15, -5, 5,
                15, 20, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 40, 45, 50, 51])
roAl = np.array([0.0, 3.803, 5.287, 7.206, 8.480, 9.375, 9.994, 10.386, 10.6,
                 10.6, 10.6, 10.6, 10.6, 10.6, 10.6, 10.6, 10.386, 9.994,
                 9.375, 8.480, 7.206, 5.747, 4.708, 4.015, 3.621, 3.5, 3.5,
                 3.5, 3.5])


def radius_of(z, zv, rv):
    """Piecewise-linear polycone outer radius; 0 outside the z range."""
    r = np.interp(z, zv, rv, left=0.0, right=0.0)
    r[(z < zv[0]) | (z > zv[-1])] = 0.0
    return r


# ── cross sections ──────────────────────────────────────────────────────────
def sigma_al_capture(E_eV):
    """27Al(n,g), 1/v below the 5.9 keV resonance region [cm2]."""
    return SIG0_AL * np.sqrt(E0 / np.asarray(E_eV))


def he3_np_interp():
    with h5py.File(ROOT / "data/He3.h5", "r") as f:
        g = f["He3"]
        E = g["energy/294K"][...]
        ds = g["reactions/reaction_103/294K/xs"]
        xs = ds[...]
        i0 = ds.attrs.get("threshold_idx", 0)
    E, xs = np.asarray(E[i0:i0 + len(xs)]), np.asarray(xs)

    def f_sig(Eq):
        lx, ly = np.log(E), np.log(np.clip(xs, 1e-300, None))
        return np.exp(np.interp(np.log(Eq), lx, ly)) * 1e-24  # cm2
    return f_sig


# ── beam ────────────────────────────────────────────────────────────────────
def load_beam():
    h = uproot.open(ROOT / "data/fluxEAR2-Ph3_in_different_units.root")[
        "flux_n_pulse_NOisolet_100bpd"]
    fe, fv = h.axis().edges(), h.values()
    g = uproot.open(ROOT / "data/lamda2DvsEn_EAR2.root")["Lambda2D"]
    lam = g.values()                                    # (240 E, 3000 r)
    lam_le = np.log10(g.axis(0).centers())
    r_edges = g.axis(1).edges()                         # cm, 10 um bins
    # rebin radius 3000 -> 60 x 0.5 mm exactly as the gun (rebin=50)
    lam60 = lam.reshape(lam.shape[0], 60, 50).sum(axis=2)
    good = lam60.sum(axis=1) > 0
    return fe, fv, lam_le[good], lam60[good], r_edges[::50]


def window_bins(fe, fv, elo, ehi):
    """Mirror the gun: every flux bin overlapping the window enters with FULL
    weight; the representative E range is clipped to the window."""
    sel = (fe[1:] > elo) & (fe[:-1] < ehi) & (fv > 0)
    lo = np.maximum(fe[:-1][sel], elo)
    hi = np.minimum(fe[1:][sel], ehi)
    return lo, hi, fv[sel]


# ── the transport integral ──────────────────────────────────────────────────
def single_pass(elo, ehi, valve_first=True, dz_mm=0.02, nsub_r=10):
    """Straight-line attenuation along the beam for every (E-bin, r).

    Returns dict with per-neutron capture probabilities and diagnostics.
    P_cap(r,E) = sum_z dz Sig_g(E) [Al] exp(-Sig_g cumAl - Sig_He3 cumGas).
    """
    fe, fv, lam_le, lam60, r60 = load_beam()
    blo, bhi, bw = window_bins(fe, fv, elo, ehi)
    phi_pulse = bw.sum()

    zg = np.arange(-35.05, 51.05 + 1e-9, dz_mm)         # mm
    ro_g = radius_of(zg, zGas, roGas)
    ro_a = radius_of(zg, zAl, roAl)

    # sub-radii: uniform within each 0.5 mm bin (as the gun samples)
    r_sub = (r60[:-1, None] + (np.arange(nsub_r)[None, :] + 0.5)
             * np.diff(r60)[:, None] / nsub_r).ravel() * 10.0   # mm
    keep = r_sub < 11.0                                  # beyond: no Al on path
    rk = r_sub[keep]

    in_gas = rk[:, None] < ro_g[None, :]
    in_al = (~in_gas) & (rk[:, None] < ro_a[None, :])

    order = slice(None, None, -1) if valve_first else slice(None)
    gas_o, al_o = in_gas[:, order], in_al[:, order]
    dz_cm = dz_mm / 10.0
    cum_gas = np.cumsum(gas_o, axis=1) * dz_cm - gas_o * dz_cm   # exclusive
    cum_al = np.cumsum(al_o, axis=1) * dz_cm - al_o * dz_cm

    t_pre = np.zeros(rk.size)                            # Al before 1st gas [cm]
    first_gas = np.argmax(gas_o, axis=1)
    has_gas = gas_o.any(axis=1)
    t_pre[has_gas] = cum_al[np.arange(rk.size), first_gas][has_gas]
    t_pre[~has_gas] = al_o[~has_gas].sum(axis=1) * dz_cm
    t_tot = al_o.sum(axis=1) * dz_cm

    # effective 1/v mean of sqrt(E0/E) over a log-uniform bin
    lnr = np.log(bhi / blo)
    m_invv = np.where(lnr > 1e-12,
                      2.0 * (blo**-0.5 - bhi**-0.5) / lnr, blo**-0.5) * np.sqrt(E0)
    Ec = np.sqrt(blo * bhi)
    sig_he3 = he3_np_interp()(Ec)

    P_al = np.zeros(len(bw))          # per-neutron capture prob, per E bin
    P_al_unshadowed = np.zeros(len(bw))
    P_he3 = np.zeros(len(bw))
    zprof = np.zeros(zg.size)         # capture density along local z [per n/bin]

    # radial pdf row per E bin (nearest log10 E row, as the gun)
    row_idx = np.abs(lam_le[None, :] - np.log10(Ec)[:, None]).argmin(axis=1)
    pdf60 = lam60 / lam60.sum(axis=1, keepdims=True)

    for b in range(len(bw)):
        w_r = np.repeat(pdf60[row_idx[b]], nsub_r)[keep] / nsub_r
        Sg = N_AL * SIG0_AL * m_invv[b]                  # 1/cm, bin-avg
        Sh = N_HE3 * sig_he3[b]
        att = np.exp(-Sg * cum_al - Sh * cum_gas)
        cap = Sg * dz_cm * al_o * att                    # (r, z_along_path)
        P_r = cap.sum(axis=1)
        P_al[b] = w_r @ P_r
        P_al_unshadowed[b] = w_r @ (cap * (cum_gas < 1e-6)).sum(axis=1)
        P_he3[b] = w_r @ (Sh * dz_cm * gas_o * att).sum(axis=1)
        zprof += (w_r @ cap)[order] * bw[b]

    res = dict(blo=blo, bhi=bhi, bw=bw, phi_pulse=phi_pulse,
               P_al=P_al, P_al_unshadowed=P_al_unshadowed, P_he3=P_he3,
               zg=zg, zprof=zprof, r_mm=rk, t_pre_cm=t_pre, t_tot_cm=t_tot,
               m_invv=m_invv, Ec=Ec)
    # volume sanity checks (full r range, not just keep)
    zg_f = np.arange(-35.0, 51.0, 0.005)
    rg = np.arange(0.0025, 12.0, 0.005)
    gasm = rg[:, None] < radius_of(zg_f, zGas, roGas)[None, :]
    alm = (~gasm) & (rg[:, None] < radius_of(zg_f, zAl, roAl)[None, :])
    ring = 2 * np.pi * rg * 0.005 * 0.005 / 1000.0       # cm3 per cell
    res["V_al_cm3"] = float((alm * ring[:, None]).sum())
    res["V_gas_cm3"] = float((gasm * ring[:, None]).sum())
    return res


def rates(res, sub_ehi=None):
    """Per-neutron, per-pulse, per-day capture rates in a sub-window."""
    m = np.ones(len(res["bw"]), bool)
    if sub_ehi is not None:
        m = res["blo"] < sub_ehi
    per_pulse = float((res["bw"] * res["P_al"])[m].sum())
    per_n = per_pulse / res["phi_pulse"]
    return per_n, per_pulse, per_pulse * PULSES_DAY


def main():
    print("=" * 72)
    print("Volume checks + geometry")
    res_v = single_pass(1e-3, 2.0, valve_first=True)
    print(f"  V_Al  = {res_v['V_al_cm3']:.3f} cm3 (target 4.907)"
          f"   V_gas = {res_v['V_gas_cm3']:.3f} cm3 (target 17.003)")

    res_n = single_pass(1e-3, 2.0, valve_first=False)

    print("\nSim window [1 meV, 2 eV], phi = %.4e n/pulse" % res_v["phi_pulse"])
    for tag, r in [("VALVE-first (as built)", res_v), ("nose-first", res_n)]:
        pn, pp, pd = rates(r)
        pn_uns = float((r["bw"] * r["P_al_unshadowed"]).sum()) / r["phi_pulse"]
        phe3 = float((r["bw"] * r["P_he3"]).sum()) / r["phi_pulse"]
        print(f"  {tag:24s} P_Al={pn:.4e}/n = {pp:9.0f}/pulse = {pd:.3e}/day"
              f"   (unshadowed part {pn_uns:.4e}; He3 abs {phe3:.3f})")

    print("\nGeant4 target: 7.806e-3 /n = 33508/pulse (He3Cap_Al nCapture)")

    for ehi, name in [(1.0, "<1 eV"), (0.5, "<0.5 eV")]:
        pn, pp, pd = rates(res_v, ehi)
        print(f"  valve-first {name:7s}: {pn:.4e}/n  {pp:9.0f}/pulse  {pd:.3e}/day")

    # extended window (note-level, resonances above ~2 keV not modelled)
    res_k = single_pass(1e-3, 1e3, valve_first=True)
    pn, pp, pd = rates(res_k)
    print(f"  valve-first <1 keV : {pn:.4e}/n  {pp:9.0f}/pulse  {pd:.3e}/day"
          f"   (phi={res_k['phi_pulse']:.4e}/pulse)")

    # ── thin-disk variants (the suspect calc), <1 eV window, sigma_th ──────
    fe, fv, lam_le, lam60, r60 = load_beam()
    _, _, bw1 = window_bins(fe, fv, 1e-3, 1.0)
    phi1 = bw1.sum()
    n_atoms = 13.24 * NA / A_AL
    A_face = np.pi * 1.06**2                             # cm2, Al outer radius
    # flux fraction inside the capsule face (thermal profile row at 25 meV)
    row = np.abs(lam_le - np.log10(0.0253)).argmin()
    pdf = lam60[row] / lam60[row].sum()
    r_cm = 0.5 * (r60[:-1] + r60[1:])
    f_hit = pdf[r_cm < 1.06].sum()
    v1 = phi1 * f_hit * (n_atoms / A_face) * SIG0_AL          # mass smeared
    v2 = phi1 * f_hit * N_AL * 0.06 * SIG0_AL                 # 0.6 mm wall
    v3 = phi1 * f_hit * N_AL * 0.55 * SIG0_AL                 # 5.5 mm nose
    m1 = float(np.mean(res_v["m_invv"][res_v["blo"] < 1.0]))  # unweighted feel
    print("\nThin-disk variants (<1 eV, sigma_th, f_hit=%.3f):" % f_hit)
    print(f"  13.24 g smeared over face (t=1.39 cm): {v1:9.0f}/pulse")
    print(f"  0.6 mm wall slab                     : {v2:9.0f}/pulse")
    print(f"  5.5 mm nose slab                     : {v3:9.0f}/pulse")
    print(f"  (flux-avg 1/v boost <sigma>/sigma_th over <1eV: "
          f"{float((bw1 * res_v['m_invv'][res_v['blo']<1.0]).sum()/phi1):.3f})")

    np.savez(OUT / "analytic_results.npz",
             **{k: v for k, v in res_v.items()},
             zprof_nose=res_n["zprof"], P_al_nose=res_n["P_al"],
             blo_k=res_k["blo"], bhi_k=res_k["bhi"], bw_k=res_k["bw"],
             P_al_k=res_k["P_al"], phi_pulse_k=res_k["phi_pulse"],
             thin_disk=[v1, v2, v3], f_hit=f_hit, phi1eV=phi1)
    json.dump({"valve_first_per_n": rates(res_v)[0],
               "nose_first_per_n": rates(res_n)[0],
               "geant4_per_n": 7.806e-3}, open(OUT / "summary.json", "w"),
              indent=1)
    print("\nsaved ->", OUT / "analytic_results.npz")


if __name__ == "__main__":
    main()
