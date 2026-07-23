#!/usr/bin/env python3
"""Figures for the Al(n,gamma) analytic-vs-Geant4 cross-check
(docs/al_gamma_yield_check/). Reads analytic_results.npz (from
al_capture_crosscheck.py) and caps_20files.root (EventTree capture count,
20 files x 1e7 neutrons of neutrons_thermal_trig_2cm).

Beam coordinate: s = y_world; capsule local z maps to y_world = -z
(placement rotateX(-90): valve local z=+51 sits at y=-51, facing the gun).
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot

ROOT = Path(__file__).resolve().parent.parent.parent   # repo root
OUT = ROOT / "docs/al_gamma_yield_check"
FIG = OUT / "figs"
FIG.mkdir(parents=True, exist_ok=True)

# palette (validated, light surface)
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#6b6a66"
ALGRAY, GASBLUE, CFRPDK = "#b9b7b2", "#d7e8f9", "#4a4945"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150})

NSIM = 2e8
PHI = 4.2924e6            # n/pulse in [1 meV, 2 eV]
PULSES_DAY = 1.929e4

d = np.load(OUT / "analytic_results.npz")
f = uproot.open(OUT / "caps_20files.root")

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
roCFRP = np.array([0.0, 4.703, 6.187, 8.106, 9.380, 10.275, 10.894, 11.286,
                   11.5, 11.5, 11.5, 11.5, 11.5, 11.5, 11.5, 11.5, 11.286,
                   10.894, 10.275, 9.380, 8.106, 6.647, 5.608, 4.915, 4.521,
                   4.4, 4.4, 4.4, 4.4])


def beam_x(z_local):
    return -np.asarray(z_local)          # y_world


# ── fig 1: geometry + Al thickness profile ──────────────────────────────────
def fig_geometry():
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.2, 3.6),
                                  gridspec_kw={"width_ratios": [1.5, 1]})
    for rv, col, lab in [(roCFRP, CFRPDK, "CFRP wrap (0.9 mm)"),
                         (roAl, ALGRAY, "Al vessel (13.24 g)"),
                         (roGas, GASBLUE, "He-3, 500 bar")]:
        x = beam_x(zAl if len(rv) == len(zAl) else zGas)
        z = zAl if len(rv) == len(zAl) else zGas
        ax.fill_between(beam_x(z), rv, -rv, color=col, lw=0, label=lab)
    ax.annotate("", xy=(-51, 14.5), xytext=(-62, 14.5),
                arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2.2))
    ax.text(-62, 16.2, "beam", color=BLUE, fontsize=11)
    ax.text(-46, -16.5, "valve + neck\n(faces beam!)", ha="center", fontsize=9)
    ax.text(30, -16.5, "nose\n(downstream)", ha="center", fontsize=9)
    ax.set_xlabel("beam coordinate  $y_{\\rm world}=-z_{\\rm capsule}$  [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_xlim(-65, 45); ax.set_ylim(-19, 19)
    ax.legend(loc="upper right", fontsize=8.5, frameon=False)

    r = d["r_mm"]
    ax2.step(r, d["t_pre_cm"] * 10, where="mid", color=BLUE, lw=2,
             label="Al seen before He-3 (as built)")
    ax2.step(r, d["t_tot_cm"] * 10, where="mid", color=GRAY, lw=1.4, ls="--",
             label="total Al on chord")
    ax2.axvspan(0, 3.5, color=BLUE, alpha=0.07, lw=0)
    ax2.text(1.7, 20.5, "beam core\n($r_{50}$=2.6 mm)", fontsize=8.5,
             ha="center", color=BLUE)
    ax2.set_xlabel("neutron radius r [mm]")
    ax2.set_ylabel("Al path [mm]")
    ax2.set_xlim(0, 11.5); ax2.set_ylim(0, 47)
    ax2.legend(fontsize=8.5, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "geometry.pdf"); plt.close(fig)


# ── fig 2: buildup of estimates ─────────────────────────────────────────────
def fig_buildup(g4_sub1eV):
    v_mass, v_wall, v_nose = d["thin_disk"]
    single_1ev = float((d["bw"] * d["P_al"])[d["blo"] < 1.0].sum())
    nose_1ev = float((d["bw"] * d["P_al_nose"])[d["blo"] < 1.0].sum())
    rows = [
        ("thin disk, 0.6 mm wall,  $\\sigma_{\\rm th}$", v_wall, GRAY),
        ("thin disk, 5.5 mm nose,  $\\sigma_{\\rm th}$", v_nose, GRAY),
        ("thin disk, 13.24 g$\\,/\\,$face,  $\\sigma_{\\rm th}$", v_mass, ORANGE),
        ("single-pass, nose into beam", nose_1ev, GRAY),
        ("single-pass, valve into beam (as built)", single_1ev, BLUE),
        ("Geant4 (full transport)", g4_sub1eV, AQUA),
    ]
    fig, ax = plt.subplots(figsize=(8.6, 3.4))
    ypos = np.arange(len(rows))[::-1]
    for y, (lab, v, c) in zip(ypos, rows):
        ax.barh(y, v, color=c, height=0.62)
        ax.text(v * 1.07, y, f"{v:,.0f}  ({v / g4_sub1eV:.2f}$\\times$)",
                va="center", fontsize=9)
        ax.text(2.0e3 * 0.92, y, lab, va="center", ha="right", fontsize=9.5)
    ax.set_xscale("log"); ax.set_xlim(2e3, 3.4e5)
    ax.set_yticks([])
    ax.set_xlabel("$^{27}$Al(n,$\\gamma$) captures / pulse   (E < 1 eV)")
    ax.spines["left"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "buildup.pdf"); plt.close(fig)
    return single_1ev, nose_1ev


# ── fig 3: capture axial profile, analytic vs Geant4 ────────────────────────
def fig_axial():
    hy = f["h_capy_al"]
    e = hy.axis().edges(); v = hy.values(flow=False)
    scale = PHI / NSIM / np.diff(e)          # per pulse per mm
    zg = d["zg"]; dzmm = zg[1] - zg[0]
    fig, ax = plt.subplots(figsize=(8.6, 3.6))
    ax.stairs(v * scale, e, fill=True, color=AQUA, alpha=0.35, lw=0)
    ax.stairs(v * scale, e, color=AQUA, lw=1.6, label="Geant4 (full transport)")
    ax.plot(beam_x(zg), d["zprof"] / dzmm, color=BLUE, lw=1.8,
            label="analytic single-pass (valve first)")
    # counterfactual flipped mounting: nose face upstream; its own beam frame
    # (beam still from the left, nose face at z=-35)
    ax.plot(zg, d["zprof_nose"] / dzmm, color=GRAY, lw=1.4,
            ls="--", label="analytic if nose faced beam (flipped)")
    ax.set_yscale("log"); ax.set_ylim(0.5, 9e3)
    ax.set_xlim(-58, 42)
    for x0, lab in [(-45.5, "valve"), (-30, "shoulder"), (0, "barrel"),
                    (32, "nose")]:
        ax.text(x0, 4.5e3, lab, ha="center", fontsize=9, color=GRAY)
    ax.axvline(-40, color=GRAY, lw=0.6, ls=":")
    ax.axvline(-21, color=GRAY, lw=0.6, ls=":")
    ax.axvline(21, color=GRAY, lw=0.6, ls=":")
    ax.set_xlabel("beam coordinate $y_{\\rm world}$ [mm]  (beam from the left)")
    ax.set_ylabel("Al captures / pulse / mm")
    ax.legend(fontsize=9, frameon=True, framealpha=0.95, edgecolor="none", loc="center right")
    fig.tight_layout()
    fig.savefig(FIG / "axial_profile.pdf"); plt.close(fig)


# ── fig 4: 2D capture map with capsule outline ──────────────────────────────
def fig_map():
    h2 = f["h_zr_al"]
    ye = h2.axis(0).edges(); re = h2.axis(1).edges()
    v = h2.values()
    fig, ax = plt.subplots(figsize=(8.6, 3.2))
    pc = ax.pcolormesh(ye, re, (v.T * PHI / NSIM), cmap="Blues",
                       norm=matplotlib.colors.LogNorm(vmin=1e-2))
    for z, rv, c in [(zAl, roAl, "k"), (zGas, roGas, BLUE)]:
        ax.plot(beam_x(z), rv, color=c, lw=1.0)
    ax.annotate("", xy=(-53, 13.2), xytext=(-59, 13.2),
                arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2))
    ax.text(-58.5, 14.0, "beam", color=BLUE, fontsize=9)
    ax.set_xlim(-60, 40); ax.set_ylim(0, 16)
    ax.set_xlabel("beam coordinate $y_{\\rm world}$ [mm]")
    ax.set_ylabel("r [mm]")
    cb = fig.colorbar(pc, ax=ax, pad=0.01)
    cb.set_label("captures / pulse / bin", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "capture_map.pdf"); plt.close(fig)


# ── fig 5: energy spectrum of Al captures ───────────────────────────────────
def fig_espec():
    hE = f["h_E_al"]
    e = hE.axis().edges(); v = hE.values()
    scale = PHI / NSIM
    # analytic into same bins
    lec = np.log10(d["Ec"])
    an = np.zeros(len(v))
    for le, w in zip(lec, d["bw"] * d["P_al"]):
        k = np.searchsorted(e, le) - 1
        if 0 <= k < len(an):
            an[k] += w
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8.6, 4.6), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1]})
    ax.stairs(v * scale, e, fill=True, color=AQUA, alpha=0.35, lw=0)
    ax.stairs(v * scale, e, color=AQUA, lw=1.6, label="Geant4")
    ax.stairs(an, e, color=BLUE, lw=1.8, label="analytic single-pass")
    ax.set_ylabel("Al captures / pulse / bin")
    ax.set_yscale("log"); ax.set_ylim(1, 3e3)
    ax.legend(frameon=False, fontsize=9)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(an > 0, v * scale / an, np.nan)
        err = np.where(an > 0, np.sqrt(v) * scale / an, np.nan)
    cx = 0.5 * (e[1:] + e[:-1])
    axr.errorbar(cx, ratio, yerr=err, fmt="o", ms=2.5, color=GRAY, lw=0.8)
    axr.axhline(1, color=GRAY, lw=0.7, ls=":")
    axr.set_ylim(0.55, 1.35)
    axr.set_ylabel("G4 / analytic", fontsize=9)
    axr.set_xlabel("log$_{10}$(E$_n$ / eV)")
    fig.tight_layout()
    fig.savefig(FIG / "espec.pdf"); plt.close(fig)


# ── numbers for the note ────────────────────────────────────────────────────
def numbers():
    hE = f["h_E_al"]
    e = hE.axis().edges(); v = hE.values(flow=True)
    n_al = v.sum()
    sub1 = v[1:-1][0.5 * (e[1:] + e[:-1]) < 0].sum() + v[0]
    hy = f["h_capy_al"]; ye = hy.axis().edges(); yv = hy.values()
    yc = 0.5 * (ye[1:] + ye[:-1])
    g4frac = {}
    for lab, lo, hi in [("valve  z 40..51", -51.5, -40), ("shoulder z 21..40", -40, -21),
                        ("barrel |z|<21", -21, 21), ("nose z<-21", 21, 36)]:
        g4frac[lab] = yv[(yc > lo) & (yc < hi)].sum() / yv.sum()
    # analytic fractions
    zg = d["zg"]; zp = d["zprof"]
    anfrac = {}
    for lab, lo, hi in [("valve  z 40..51", 40, 51.5), ("shoulder z 21..40", 21, 40),
                        ("barrel |z|<21", -21, 21), ("nose z<-21", -36, -21)]:
        anfrac[lab] = zp[(zg > lo) & (zg < hi)].sum() / zp.sum()
    print("Geant4 Al captures (20 files): %d  -> %.4e /n" % (n_al, n_al / NSIM))
    print("  <1 eV fraction: %.4f" % (sub1 / n_al))
    print("  region fractions  (Geant4  |  analytic):")
    for k in g4frac:
        print(f"    {k:20s} {g4frac[k]:6.3f}  |  {anfrac[k]:6.3f}")
    return n_al / NSIM * PHI * (sub1 / n_al)


g4_sub1eV = numbers()
fig_geometry()
s1, n1 = fig_buildup(g4_sub1eV)
fig_axial()
fig_map()
fig_espec()
print("single-pass <1eV: %.0f   nose <1eV: %.0f   G4 <1eV: %.0f" %
      (s1, n1, g4_sub1eV))
print("figs ->", FIG)
