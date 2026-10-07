#!/usr/bin/env python3
"""scattering_geometry.py — acceptance-averaged upstream material budgets.

Ray-traces e+e- tracks from the He-3 gas through the *actual* sim geometry and
into the MM acceptance, instead of assuming normal incidence through the
barrel.  Per track it accumulates x/X0 of every material crossed, so nose/butt
exits, wall obliquity and the oblique path through the air gap / window / drift
gas are all included.

  baseline : STEP capsule profile copied from src/DetectorConstruction.cc
             (gas / Al / CFRP polycones), vertices uniform in the gas.
  upgrade  : 2 bar mylar-wrap cell from x17_facility_search/vessel_design
             (hexagonal 12 um aluminised-mylar skin, rods/seam ignored).
             Tracks leaving through an end cap are lost (cell end caps are
             thick) and are counted as acceptance loss, not scattering.

MM acceptance: four flat arms at distance D from the beam axis, active area
u = 399 mm (tangential) x v = 360 mm (along beam), centred on the origin
(pinwheel shifts ignored).  Beam axis = z.
"""
import numpy as np

# radiation lengths [cm] ------------------------------------------------------
X0 = dict(Al=8.897, CFRP=42.70 / 1.55, Mylar=28.54, Kapton=28.58, Cu=1.436,
          air=36.62 / 1.205e-3, He4=94.32 / 1.66e-4, ArIso=10900.0)


def he3_x0_cm(p_bar, T=293.15):
    rho = p_bar * 101325 * 3.016e-3 / (8.314 * T) / 1000
    return 70.7 / rho


# ---------------------------------------------------------------- capsule ----
Z_GAS = np.array([-29.5, -28, -26, -24, -22, -20, -15, -5, 5, 15, 20, 22, 24, 26,
                  28, 30, 32, 34, 36, 38, 40, 44, 50.7])
R_GAS = np.array([0.001, 6, 8, 9.165, 9.798, 10, 10, 10, 10, 10, 10, 9.798, 9.165,
                  8, 6.299, 4.842, 3.66, 2.711, 1.967, 1.41, 1.026, 0.75, 0.75])
Z_V = np.array([-35, -34, -33, -31, -29, -27, -25, -23, -21, -20, -15, -5, 5, 15,
                20, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 40, 45, 50, 51.])
R_AL = np.array([0, 3.803, 5.287, 7.206, 8.48, 9.375, 9.994, 10.386, 10.6, 10.6,
                 10.6, 10.6, 10.6, 10.6, 10.6, 10.6, 10.386, 9.994, 9.375, 8.48,
                 7.206, 5.747, 4.708, 4.015, 3.621, 3.5, 3.5, 3.5, 3.5])
R_CF = np.array([0, 4.703, 6.187, 8.106, 9.38, 10.275, 10.894, 11.286, 11.5, 11.5,
                 11.5, 11.5, 11.5, 11.5, 11.5, 11.5, 11.286, 10.894, 10.275, 9.38,
                 8.106, 6.647, 5.608, 4.915, 4.521, 4.4, 4.4, 4.4, 4.4])


def _prof(z, zt, rt):
    return np.where((z >= zt[0]) & (z <= zt[-1]), np.interp(z, zt, rt), -1.0)


def _directions(n, rng):
    c = rng.uniform(-1, 1, n)
    s = np.sqrt(1 - c * c)
    ph = rng.uniform(0, 2 * np.pi, n)
    return np.stack([s * np.cos(ph), s * np.sin(ph), c], 1)


def mm_acceptance(v, d, D=250.0, u_hw=199.5, v_hw=180.0):
    """bool mask + t_plane + |cos| to the plane actually hit (mm)."""
    n = len(v)
    ok = np.zeros(n, bool)
    t_hit = np.zeros(n)
    cosp = np.ones(n)
    for nx, ny in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        dn = d[:, 0] * nx + d[:, 1] * ny
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (D - (v[:, 0] * nx + v[:, 1] * ny)) / dn
        p = v + t[:, None] * d
        u = p[:, 0] * (-ny) + p[:, 1] * nx
        hit = (dn > 1e-6) & (np.abs(p[:, 2]) <= v_hw) & (np.abs(u) <= u_hw) & ~ok
        t_hit[hit] = t[hit]
        cosp[hit] = dn[hit]
        ok |= hit
    return ok, t_hit, cosp


def _sample_gas_vertices(n, rng):
    out = []
    while sum(len(o) for o in out) < n:
        q = np.stack([rng.uniform(-10, 10, 4 * n), rng.uniform(-10, 10, 4 * n),
                      rng.uniform(-29.5, 50.7, 4 * n)], 1)
        ok = np.hypot(q[:, 0], q[:, 1]) < np.interp(q[:, 2], Z_GAS, R_GAS)
        out.append(q[ok])
    return np.concatenate(out)[:n]


def baseline(n=300_000, seed=1, D=250.0, ds=0.1):
    rng = np.random.default_rng(seed)
    v = _sample_gas_vertices(3 * n, rng)
    d = _directions(len(v), rng)
    ok, tpl, cosp = mm_acceptance(v, d, D)
    v, d, tpl, cosp = v[ok][:n], d[ok][:n], tpl[ok][:n], cosp[ok][:n]
    n = len(v)
    L = dict(gas=np.zeros(n), Al=np.zeros(n), CFRP=np.zeros(n))
    pos = v.copy()
    alive = np.ones(n, bool)
    t_exit = np.zeros(n)
    for i in range(1, 1200):
        pos = v + (i * ds) * d
        r = np.hypot(pos[:, 0], pos[:, 1]); z = pos[:, 2]
        in_gas = r < _prof(z, Z_GAS, R_GAS)
        in_al = (~in_gas) & (r < _prof(z, Z_V, R_AL))
        in_cf = (~in_gas) & (~in_al) & (r < _prof(z, Z_V, R_CF))
        inside = in_gas | in_al | in_cf
        for k, m in (("gas", in_gas), ("Al", in_al), ("CFRP", in_cf)):
            L[k] += np.where(alive & m, ds, 0.0)
        newly = alive & ~inside
        t_exit[newly] = i * ds
        alive &= inside
        if not alive.any():
            break
    t_exit[alive] = 1200 * ds
    z_exit = v[:, 2] + t_exit * d[:, 2]
    nose_butt = np.abs(z_exit) > 20.0          # beyond the barrel: dome / neck
    cm = 0.1
    x = dict(
        he3=L["gas"] * cm / he3_x0_cm(500.0),
        Al=L["Al"] * cm / X0["Al"],
        CFRP=L["CFRP"] * cm / X0["CFRP"],
        air=np.maximum(tpl - t_exit, 0) * cm / X0["air"],
        win_Mylar=0.004 / cosp / X0["Mylar"],
        win_Kapton=0.005 / cosp / X0["Kapton"],
        win_Cu=0.0009 / cosp / X0["Cu"],
        drift=3.0 / cosp / X0["ArIso"],
    )
    xcap = x["he3"] + x["Al"] + x["CFRP"]
    zb = np.arange(-36, 54, 2.0)
    ib = np.clip(np.digitize(z_exit, zb) - 1, 0, len(zb) - 2)
    nb = np.bincount(ib, minlength=len(zb) - 1)
    xs = np.bincount(ib, weights=xcap, minlength=len(zb) - 1)
    diag = dict(zbins=zb.tolist(), n_by_bin=nb.tolist(),
                x_by_bin=(xs / np.maximum(nb, 1)).tolist(),
                n=n, nose_butt_frac=float(nose_butt.mean()),
                al_mm_barrel=float(L["Al"][~nose_butt].mean()),
                al_mm_nosebutt=float(L["Al"][nose_butt].mean()),
                x_tot_barrel=float(sum(x[k][~nose_butt] for k in ("he3", "Al", "CFRP")).mean()),
                x_tot_nosebutt=float(sum(x[k][nose_butt] for k in ("he3", "Al", "CFRP")).mean()))
    return x, diag


def upgrade(n=300_000, seed=2, D=250.0, p_bar=2.0, t_mylar=12e-4, t_al=2 * 40e-7,
            n_extra_foils=0, he4_gap=False, win_foils=2, t_win=12e-4,
            apothem=34.775, half_len=150.0, r_beam=30.0):
    """2 bar mylar-wrap cell.  Lengths in mm (cm for foil thicknesses)."""
    rng = np.random.default_rng(seed)
    N = 4 * n
    v = np.stack([rng.uniform(-r_beam, r_beam, N), rng.uniform(-r_beam, r_beam, N),
                  rng.uniform(-half_len, half_len, N)], 1)
    v = v[np.hypot(v[:, 0], v[:, 1]) <= r_beam]
    d = _directions(len(v), rng)
    ok, tpl, cosp = mm_acceptance(v, d, D)
    v, d, tpl, cosp = v[ok], d[ok], tpl[ok], cosp[ok]
    # hexagonal skin, facet normals between the rods (rods at pi/2 + k*60 deg)
    phis = np.pi / 2 + np.pi / 6 + np.arange(6) * np.pi / 3
    t_fac = np.full(len(v), np.inf); cosn = np.ones(len(v))
    for ph in phis:
        nx, ny = np.cos(ph), np.sin(ph)
        dn = d[:, 0] * nx + d[:, 1] * ny
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (apothem - (v[:, 0] * nx + v[:, 1] * ny)) / dn
        t = np.where(dn > 1e-9, t, np.inf)
        better = t < t_fac
        t_fac[better] = t[better]; cosn[better] = dn[better]
    with np.errstate(divide="ignore", invalid="ignore"):
        t_end = np.where(np.abs(d[:, 2]) > 1e-9,
                         (np.sign(d[:, 2]) * half_len - v[:, 2]) / d[:, 2], np.inf)
    lost = t_end < t_fac
    keep = ~lost
    frac_lost = float(lost.mean())
    v, d, tpl, cosp, t_fac, cosn = (a[keep][:n] for a in (v, d, tpl, cosp, t_fac, cosn))
    cm = 0.1
    xh = he3_x0_cm(p_bar)
    x = dict(he3=t_fac * cm / xh,
             Mylar_cell=t_mylar / cosn / X0["Mylar"],
             Al_cell=t_al / cosn / X0["Al"],
             win_Mylar=win_foils * t_win / cosp / X0["Mylar"],
             win_Al=win_foils * (2 * 40e-7) / cosp / X0["Al"],
             drift=3.0 / cosp / X0["ArIso"])
    gap = np.maximum(tpl - t_fac, 0) * cm
    if he4_gap:
        x["he4"] = gap / X0["He4"]
        x["bag_Mylar"] = n_extra_foils * t_mylar / cosp / X0["Mylar"]
        x["bag_Al"] = n_extra_foils * (2 * 40e-7) / cosp / X0["Al"]
    else:
        x["air"] = gap / X0["air"]
    diag = dict(n=len(v), frac_lost_endcap=frac_lost,
                path_factor_cell=float((1 / cosn).mean()),
                path_factor_plane=float((1 / cosp).mean()))
    return x, diag


def mean(x):
    return {k: float(np.mean(v)) for k, v in x.items()}


if __name__ == "__main__":
    xb, db = baseline()
    print("baseline diag", db)
    mb = mean(xb)
    print({k: round(v * 100, 4) for k, v in mb.items()}, "tot upstream %",
          round(sum(v for k, v in mb.items() if k != 'drift') * 100, 3))
    xu, du = upgrade()
    print("upgrade diag", du)
    mu = mean(xu)
    print({k: round(v * 100, 4) for k, v in mu.items()}, "tot upstream %",
          round(sum(v for k, v in mu.items() if k != 'drift') * 100, 3))
