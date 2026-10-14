#!/usr/bin/env python3
"""make_scattering_upgrade_deck.py — slide note: how much does cutting the
material budget sharpen the X17 opening-angle peak?

Everything (budgets, scattering curves, peak shapes) is recomputed from
scattering_geometry.py + make_scattering_upgrade_figs.py + the Geant4 response,
so rerunning after a change updates numbers, charts and claims together.

Usage:  python scripts/make_scattering_upgrade_deck.py [out.html]
"""
import sys
import os
import json
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.expanduser(os.environ.get(
    'SLIDEDOC_DIR', '~/PycharmProjects/dylan-cern-site/scripts')))
import slidedoc as sd                       # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_scattering_upgrade_figs as F    # noqa: E402
import scattering_geometry as geo           # noqa: E402
import make_angular_resolution_figs as ar   # noqa: E402

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / "docs/scattering_upgrade/scattering-upgrade.html"

# ── numbers ──────────────────────────────────────────────────────────────────
d = json.load(open(ar.RESP))
med = F.med_sim(d)
M, DIAG = F.build_cases()
mb, mu, ma = M["base"], M["up"], M["app"]
xb, xu, xa = F.xtot(mb), F.xtot(mu), F.xtot(ma)
xb_up = xb - mb["drift"]
spec3 = {"base": ("Baseline", sd.RED, "-"), "up": ("Upgrade", sd.BLUE, "-"),
         "app": ("Upgrade + ⁴He bag", sd.ORANGE, "-")}
m7 = dict(mu); m7["Al_cell"] = mu["Al_cell"] * (7e-4 / 8e-6)    # 7 um rolled-Al laminate
M7 = dict(M); M7["up7"] = m7
cen, ht, RES = F.peak_curves(M7, d, ["base", "up", "app", "up7"])
psi = lambda m, k: med(k) * F.highland(k, F.xtot(m)) / F.highland(k, xb_up)

# colours per component, same hues as the matplotlib figures
COL = {k: v[2] for k, v in F.COMP.items()}
LAB = {k: v[1] for k, v in F.COMP.items()}
FAM = {k: v[0] for k, v in F.COMP.items()}
FAMCOL = {"³He gas": "#4a90d9", "Capsule / cell wall": "#d62728", "Gap to MM": "#2ca02c",
          "MM window": "#7a52b8", "Drift gas": "#6f7885"}
fam_share = lambda m, f: sum(v for k, v in m.items() if FAM[k] == f) / F.xtot(m)

D = sd.Deck("Scattering budget upgrade",
            "How much a 2 bar mylar 3He cell and a mylar MM window sharpen the X17 opening-angle peak")


# ── helpers ──────────────────────────────────────────────────────────────────
def pctf(v, dd=3):
    return f"{v * 100:.{dd}f}%"


def seg_tip(k, v, m):
    return (f"{LAB[k]}\nx/X₀ = {pctf(v, 4)}\n{v / F.xtot(m) * 100:.1f}% of this stack "
            f"(acceptance-averaged, incl. drift gas)")


def spread(ys, gap):
    ys = list(ys)
    for _ in range(300):
        for i in range(1, len(ys)):
            if ys[i] - ys[i - 1] < gap:
                m_ = .5 * (ys[i] + ys[i - 1]); ys[i - 1] = m_ - gap / 2; ys[i] = m_ + gap / 2
    return ys


def stack_svg_parts(m, xc, w, yb, scale, tipm=None, lab_x=None, lab_gap=26, labels=True,
                    size=21, fam_bracket_x=None):
    """Vertical stacked bar centred at xc, base at pixel yb, scale px per unit x/X0.
    Returns svg string + list of (family, y_top, y_bot) spans."""
    out, y, segs = [], yb, []
    for k, v in F.stack(m):
        h = v * scale
        out.append(f'<rect x="{xc - w / 2:.1f}" y="{y - h:.1f}" width="{w}" height="{max(h, .6):.2f}" '
                   f'fill="{COL[k]}" stroke="#fff" stroke-width=".8"{sd.tipattr(seg_tip(k, v, tipm or m))}/>')
        segs.append((k, v, y - h / 2, y - h, y))
        y -= h
    spans = {}
    for k, v, ymid, ytop, ybot in segs:
        f_ = FAM[k]
        lo, hi = spans.get(f_, (ytop, ybot))
        spans[f_] = (min(lo, ytop), max(hi, ybot))
    if labels:
        ys = spread([s[2] for s in segs][::-1], lab_gap)       # ascending y (top first)
        ys[-1] = min(ys[-1], yb - 8)
        for i in range(len(ys) - 2, -1, -1):
            ys[i] = min(ys[i], ys[i + 1] - lab_gap)
        ys = ys[::-1]
        for (k, v, ymid, _t, _b), yl in zip(segs, ys):
            out.append(sd.line(xc + w / 2, ymid, lab_x - 8, yl, "#9aa1ad", 1.2))
            txt = f"{LAB[k]}  {pctf(v, 4 if v < 1e-4 else 3)}"
            out.append(sd.T(lab_x, yl + 7, txt, size, sd.INK, "start", tip=seg_tip(k, v, tipm or m)))
    return "".join(out), spans, segs


def fam_brackets(spans, x, names=None, size=21):
    o = []
    items = sorted(spans.items(), key=lambda kv: (kv[1][0] + kv[1][1]) / 2)
    ys = spread([(lo + hi) / 2 for _, (lo, hi) in items], size + 10)
    ys[-1] = min(ys[-1], max(hi for _, (lo, hi) in items) - 4)
    for i in range(len(ys) - 2, -1, -1):
        ys[i] = min(ys[i], ys[i + 1] - (size + 10))
    for (f_, (lo, hi)), yl in zip(items, ys):
        o.append(f'<path d="M{x - 6},{lo:.1f} L{x},{lo:.1f} L{x},{hi:.1f} L{x - 6},{hi:.1f}" fill="none" '
                 f'stroke="{FAMCOL[f_]}" stroke-width="4"/>')
        o.append(sd.T(x + 12, yl + 7, f_, size, FAMCOL[f_], "start", 600))
    return "".join(o)


# ── slide 1: cover ───────────────────────────────────────────────────────────
r = RES
def cover_plot(w=1100, h=680):
    """Black card: previous vs upgraded X17 opening-angle spectrum, big type."""
    ml, mr, mt, mb = 56, 40, 96, 118
    pw_, ph_ = w - ml - mr, h - mt - mb
    ymax = .06
    X = lambda t: ml + (t - 60) / 120 * pw_
    Y = lambda v: mt + ph_ * (1 - v / ymax)
    COL_OLD, COL_NEW = "#ff7b7b", "#5fd0ff"
    o = [f'<rect width="{w}" height="{h}" rx="22" fill="{sd.DARK}"/>',
         sd.T(ml, 58, "Expected X17 opening-angle spectrum", 34, "#f2f4f7", "start", 600)]
    for t in range(60, 181, 20):
        o.append(sd.line(X(t), mt + ph_, X(t), mt + ph_ + 10, "#8b93a1", 2))
        o.append(sd.T(X(t), mt + ph_ + 46, str(t), 30, "#c9d0db"))
    o.append(sd.line(ml, mt + ph_, ml + pw_, mt + ph_, "#8b93a1", 2))
    o.append(sd.T(ml + pw_ / 2, h - 16, "e⁺e⁻ opening angle [deg]", 32, "#e6e9ee"))
    o.append(sd.line(X(109.4), mt + 10, X(109.4), mt + ph_, "#6b7280", 2, "8 7"))
    o.append(sd.T(X(109.4) - 10, mt + 36, "109° kinematic edge", 26, "#9aa3b2", "end"))
    for key, colr, nm in (("base", COL_OLD, "Previous (n_TOF)"), ("up", COL_NEW, "Upgraded")):
        xs = [c for c in cen if 60 <= c <= 180]
        sel = [i for i, c in enumerate(cen) if 60 <= c <= 180]
        ys = [min(RES[key]["h"][i], ymax) for i in sel]
        o.append(sd.poly([X(a_) for a_ in xs], [Y(b_) for b_ in ys], colr, 6.5,
                         tip=f"{nm}: σ68 = {RES[key]['s68']:.1f}°, peak ×{RES[key]['rel']:.2f} of ideal, FWHM {RES[key]['fwhm']:.0f}°"))
        for a_, b_ in list(zip(xs, ys))[::3]:
            o.append(f'<circle class="hit" cx="{X(a_):.1f}" cy="{Y(b_):.1f}" r="14" fill="transparent"'
                     f'{sd.tipattr(f"{nm}" + chr(10) + f"θ = {a_:.0f}°: {b_:.4f} per degree")}/>')
    # direct labels
    o.append(sd.T(X(131), Y(.0525), "Upgraded", 44, COL_NEW, "start", 700))
    o.append(sd.T(X(131), Y(.0525) + 40, f"σ68 {RES['up']['s68']:.1f}° scattering", 30, COL_NEW, "start"))
    o.append(sd.line(X(128), Y(.0525) - 12, X(114), Y(.0505), COL_NEW, 2.5))
    o.append(sd.T(X(131), Y(.033), "Previous (n_TOF)", 44, COL_OLD, "start", 700))
    o.append(sd.T(X(131), Y(.033) + 40, f"σ68 {RES['base']['s68']:.1f}° scattering", 30, COL_OLD, "start"))
    o.append(sd.line(X(128), Y(.033) - 12, X(119), Y(.0215), COL_OLD, 2.5))
    return sd.svg(w, h, "".join(o), "previous vs upgraded X17 opening-angle spectrum")


cover = (sd.kicker("X17 pair spectrometer · multiple scattering budget") +
         '<h1 style="font-size:80px;font-weight:600;line-height:1.08;letter-spacing:-2px;max-width:1664px">'
         f'Reducing MX17 Material Budget<br>→ {r["base"]["s68"]:.0f}° to {r["up"]["s68"]:.0f}° scattering</h1>' +
         '<div style="position:absolute;left:128px;top:420px;width:540px;display:flex;flex-direction:column;gap:56px">' +
         sd.bignum(f"{xb / xu:.0f}×", "less material", sd.DBLUE,
                   f"{xb * 100:.2f}% → {xu * 100:.2f}% X₀, drift gas included", size=120, tip=
                   "Acceptance-averaged x/X₀ over tracks that reach a MM arm, including the 30 mm drift gas.") +
         sd.bignum(f'{r["base"]["s68"]:.0f}° → {r["up"]["s68"]:.1f}°', "opening-angle σ68", sd.DRED,
                   f'median per-track scattering at 5 MeV: {psi(mb, 5):.0f}° → {psi(mu, 5):.1f}°', size=96, tip=
                   "σ68 = half the 16–84 percentile span of θ_reco − θ_truth, X17 at-rest pairs, direction at the first MM hit.") +
         '</div>' +
         f'<div style="position:absolute;left:740px;top:384px">{cover_plot()}</div>')
D.slide("cover", cover, dark=True, short="Answer",
        notes="Baseline = the STEP capsule used in pairs_v2 (500 bar ³He, 0.6 mm Al + 0.9 mm CFRP) with the standard "
              "MM window (40 µm Mylar, 50 µm Kapton, 9 µm Cu). Upgrade = 2 bar ³He in the 12 µm aluminised-mylar rod-cage "
              "cell of x17_facility_search/vessel_design, and a drift window of two aluminised mylar foils (outer for gas "
              "containment, inner for HV; 12 µm each assumed). Rough estimate: no acceptance, energy-loss or range changes.")

# ── slide 2: setup ───────────────────────────────────────────────────────────
base_flow = sd.flow([
    dict(label="³He gas", sub="500 bar, r = 10 mm bore", color="#4a90d9",
         tip="Capsule bore: r = 10 mm, ±20 mm barrel + hemispherical nose and conical neck; ρ = 62.7 mg/cm³."),
    dict(label="Al capsule", sub="0.6 mm barrel, 5.5 mm at the nose tip", color="#7f1d1d",
         tip="STEP-derived profile in DetectorConstruction.cc: Al outer radius 10.6 mm; solid Al neck r = 3.5 mm."),
    dict(label="CFRP wrap", sub="0.9 mm", color="#d62728", tip="CFRP, ρ = 1.55 g/cm³, X₀ = 27.6 cm."),
    dict(label="Air gap", sub="≈ 24 cm", color="#2ca02c", tip="Capsule surface to MM front face (250 mm from the origin in pairs_v2)."),
    dict(label="MM window + cathode", sub="Mylar 40 µm, Kapton 50 µm, Cu 9 µm", color="#7a52b8",
         tip="Drift window and cathode foil stack in front of the 30 mm drift gap."),
    dict(label="Drift gas", sub="30 mm ArIso", color="#6f7885", tip="Ar/iC₄H₁₀; X₀ ≈ 109 m. Included in every figure here."),
], size=24)
up_flow = sd.flow([
    dict(label="³He gas", sub="2 bar, beam envelope 60 mm", color="#4a90d9",
         tip="Gas path to the skin is the vertex-dependent distance inside the hexagonal cell (apothem 34.8 mm)."),
    dict(label="Mylar skin", sub="12 µm + 2×40 nm Al", color="#d62728", fill="#fdf0f0",
         tip="① Replacement 1. No pressure shell: a double-aluminised mylar sheet wrapped on a six-rod carbon cage (x17_facility_search/vessel_design). Rods and seam cover ≈ 9% of azimuth and are ignored here."),
    dict(label="Air gap", sub="≈ 24 cm (unchanged)", color="#2ca02c", tip="Unchanged in the main comparison; replaced by ⁴He in the appendix."),
    dict(label="Two mylar foils", sub="12 µm each (assumed)", color="#7a52b8", fill="#f4effa",
         tip="② Replacement 2. Outer foil contains the gas, inner foil holds the drift HV. Thickness is an assumption."),
    dict(label="Drift gas", sub="30 mm ArIso (unchanged)", color="#6f7885", tip="Unchanged."),
], size=24)
D.slide("setup",
        sd.title("Two replacements, same track path: the capsule wall and the MM window go",
                 "What an e⁺/e⁻ crosses between the ³He vertex and the first MM strips") +
        sd.p("Baseline (pairs_v2 as simulated)", 28, sd.RED, 600) + base_flow +
        sd.p("Upgrade", 28, sd.BLUE, 600) + up_flow +
        sd.callout("Hover the dotted terms, the boxes and every bar or curve for values and provenance. "
                   "Budgets are ray-traced over the actual MM acceptance (four flat arms, 399 × 360 mm), "
                   "not taken at normal incidence through the barrel.", sd.GOLD, 24),
        short="Setup",
        notes="Radiation lengths used (cm): Al 8.90, CFRP 27.6, Mylar 28.5, Kapton 28.6, Cu 1.44, air 3.04×10⁴, "
              "ArIso ≈ 1.09×10⁴, ⁴He 5.7×10⁵; ³He from the ideal-gas density (500 bar: X₀ = 1.13×10³ cm; 2 bar: 2.8×10⁵ cm). "
              "Scripts: scripts/scattering_geometry.py (ray tracing), scripts/make_scattering_upgrade_figs.py (matplotlib versions of these figures), "
              "scripts/make_scattering_upgrade_deck.py (this note).")

# ── slide 3: budget bars + exploded ──────────────────────────────────────────
W, H = 1664, 600
yb = 540
scale = 470 / xb
inner = []
s1, sp1, seg1 = stack_svg_parts(mb, 130, 150, yb, scale, lab_x=330, lab_gap=27, size=20)
inner.append(s1)
inner.append(sd.T(130, yb + 36, "Baseline", 24, sd.INK, "middle", 600))
inner.append(sd.T(130, yb - xb * scale - 14, f"{xb * 100:.2f}%", 26, sd.INK, "middle", 600))
# upgrade (true scale) + exploded
ux = 800
s2, sp2, seg2 = stack_svg_parts(mu, ux, 120, yb, scale, labels=False)
inner.append(s2)
inner.append(sd.T(ux, yb + 36, "Upgrade", 24, sd.INK, "middle", 600))
inner.append(sd.T(ux, yb - xu * scale - 14, f"{xu * 100:.2f}%", 26, sd.INK, "middle", 600))
inner.append(f'<rect x="{ux - 66}" y="{yb - xu * scale - 4:.1f}" width="132" height="{xu * scale + 8:.1f}" fill="none" stroke="#777" stroke-dasharray="5 4"/>')
ex, ew = 1010, 110
zs = 470 / xu
s3, sp3, seg3 = stack_svg_parts(mu, ex, ew, yb, zs, lab_x=1085, lab_gap=27, size=19)
inner.append(s3)
inner.append(sd.line(ux + 66, yb - xu * scale - 4, ex - ew / 2, yb - 470, "#999", 1.2, "6 5"))
inner.append(sd.line(ux + 66, yb + 4, ex - ew / 2, yb, "#999", 1.2, "6 5"))
inner.append(sd.T(ex, yb + 36, f"Upgrade, exploded ×{zs / scale:.0f}", 24, sd.INK, "middle", 600))
inner.append(fam_brackets(sp3, 1480, size=18))
# axis ticks for exploded
for t in np.arange(0, xu * 100 + 1e-9, 0.02):
    y = yb - t / 100 * zs
    inner.append(sd.line(ex - ew / 2 - 5, y, ex - ew / 2, y, "#888", 1.2))
    inner.append(sd.T(ex - ew / 2 - 10, y + 7, f"{t:.2f}%", 17, sd.MUT, "end"))
for t in (0, .5, 1.0, 1.5):
    y = yb - t / 100 * scale
    inner.append(sd.line(30, y, 38, y, "#888", 1.2))
    inner.append(sd.T(24, y + 7, f"{t:.1f}%", 17, sd.MUT, "end"))
share_wall = fam_share(mb, "Capsule / cell wall")
share_air = mu["air"] / xu
slide3 = (sd.title(
    f"The wall is {share_wall * 100:.0f}% of the baseline budget; after the upgrade {share_air * 100:.0f}% is air",
    "Acceptance-averaged x/X₀ vertex → MM strips, drift gas included; one colour family = one physical part") +
          sd.svg(W, H, "".join(inner), "stacked material budgets, baseline vs upgrade, with exploded upgrade stack"))
D.slide("budget", slide3, short="Material budget",
        foot="x/X₀ summed along each track and averaged over tracks that reach a MM arm; light shade = mylar/plastic, dark = metal.",
        notes=("Baseline components (x/X₀ %): " + ", ".join(f"{LAB[k]} {v * 100:.3f}" for k, v in F.stack(mb)) +
               ". Upgrade: " + ", ".join(f"{LAB[k]} {v * 100:.4f}" for k, v in F.stack(mu)) +
               f". The normal-incidence barrel-only budget in angular_resolution_note.md is 1.26%; the ray-traced average is "
               f"{(xb - mb['drift']) * 100:.2f}% upstream because tracks cross the wall obliquely and ~{DIAG['base']['nose_butt_frac'] * 100:.0f}% leave "
               "through the dome or neck."))

# ── slide 4: scattering vs energy by component ───────────────────────────────
ke = np.linspace(2, 15, 66)
pw, ph = 810, 560


def energy_panel(m, label, colour):
    P = sd.Plot(pw, ph, x=(2, 15), y=(0, 45), xlabel="e± kinetic energy [MeV]",
                ylabel="median scattering angle [deg]", title=f"{label}: x/X₀ = {F.xtot(m) * 100:.2f}%")
    P.xticks([(k, str(k)) for k in (2, 4, 6, 8, 10, 12, 14)])
    P.yticks([(k, str(k)) for k in (0, 10, 20, 30, 40)])
    xt = F.xtot(m)
    tot = med(ke) * F.highland(ke, xt) / F.highland(ke, xb_up)
    bot = np.zeros_like(ke)
    for k, v in F.stack(m):
        top = bot + v / xt * tot
        tip = (f"{LAB[k]}\n{v / xt * 100:.1f}% of the scattering variance\n"
               f"≈ {v / xt * psi(m, 5):.1f}° of the {psi(m, 5):.1f}° total at 5 MeV")
        P.band(list(ke), list(bot), list(top), COL[k], alpha=.95, tip=tip)
        bot = top
    P.line(list(ke), list(tot), sd.INK, w=2.5, markers=False,
           tips=[f"KE {a:.1f} MeV: total {b:.1f}°" for a, b in zip(ke, tot)])
    return P


leg4 = sd.legend([(LAB[k], COL[k], "box")
                  for k in F.ORDER if k in set(mb) | set(mu)][:14], size=19)
p4 = sd.row(energy_panel(mb, "Baseline", sd.RED).svg("scattering vs energy, baseline"),
            energy_panel(mu, "Upgrade", sd.BLUE).svg("scattering vs energy, upgrade"), gap=44)
D.slide("energy",
        sd.title(f"Per-track scattering falls ×{psi(mb, 10) / psi(mu, 10):.1f}: at 10 MeV {psi(mb, 10):.1f}° → {psi(mu, 10):.1f}°",
                 "Median space angle of the direction at the first MM hit vs e± energy; band height = share of the scattering variance (∝ x/X₀)") +
        p4 + leg4,
        short="vs energy",
        notes="Baseline total = measured Geant4 median ψ(KE) from the 'first' estimator in geant4_response.json (X17 + IPC, e⁻ and e⁺ combined), "
              "times a 1% factor for the drift gas, which that estimator does not contain. Upgrade total = the same curve times the Highland ratio "
              "θ₀(x_upgrade)/θ₀(x_baseline upstream). Bands apportion the total by each component's share of x/X₀ (scattering variances add). "
              "At X17 energies KE⁻+KE⁺ ≈ 19.6 MeV, so one leg is always ≤ 9.8 MeV; the soft leg sets the opening-angle resolution. "
              "Below 2 MeV the response table is too coarse to show.")

# ── slide 5: nose/butt exits ─────────────────────────────────────────────────
db = DIAG["base"]; du = DIAG["up"]
zb = np.array(db["zbins"]); nb = np.array(db["n_by_bin"], float); xbn = np.array(db["x_by_bin"])
PW, PH = 1000, 620
sx = lambda z: 60 + (z + 40) * (PW - 80) / 96.0       # z in [-40, 56]
cy_prof, scl = 500, 12.0                                # profile centreline (px), px per mm radius
o5 = []
# histogram of exit position above
hmax = nb.max()
for i in range(len(nb)):
    z0, z1 = zb[i], zb[i + 1]
    h = nb[i] / hmax * 190
    nose = abs(.5 * (z0 + z1)) > 20
    o5.append(f'<rect x="{sx(z0):.1f}" y="{250 - h:.1f}" width="{sx(z1) - sx(z0) - 1:.1f}" height="{h:.1f}" '
              f'fill="{sd.RED if nose else sd.GREY}" opacity=".85"{sd.tipattr(f"exit z {z0:.0f}…{z1:.0f} mm: {nb[i]:.0f} tracks ({nb[i] / nb.sum() * 100:.1f}%)" + chr(10) + f"mean capsule x/X₀ {xbn[i] * 100:.2f}%")}/>')
o5.append(sd.T(60, 56, "where accepted tracks leave the capsule (z along the beam)", 22, sd.MUT, "start"))
# x/X0 line
ok5 = [i for i in range(len(nb)) if nb[i] > 4000]
o5.append(sd.poly([sx(.5 * (zb[i] + zb[i + 1])) for i in ok5], [250 - xbn[i] / .045 * 190 for i in ok5], sd.BLUE, 3.5))
o5.append(sd.T(sx(56), 96, "blue line: capsule x/X₀ crossed", 19, sd.BLUE, "end"))
o5.append(sd.T(sx(56), 118, f"({db['x_tot_barrel'] * 100:.1f}% in the barrel, up to {max(xbn[i] for i in ok5) * 100:.1f}% at the ends)", 19, sd.BLUE, "end"))
# capsule profile
def prof(zt, rt):
    return [(sx(z), cy_prof - r * scl) for z, r in zip(zt, rt)] + [(sx(z), cy_prof + r * scl) for z, r in zip(zt[::-1], rt[::-1])]
for zt, rt, col in ((geo.Z_V, geo.R_CF, "#d62728"), (geo.Z_V, geo.R_AL, "#7f1d1d"), (geo.Z_GAS, geo.R_GAS, "#4a90d9")):
    pts = prof(list(zt), list(rt))
    o5.append(f'<polygon points="{" ".join(f"{a:.1f},{b:.1f}" for a, b in pts)}" fill="{col}"{sd.tipattr({"#d62728": "CFRP wrap 0.9 mm", "#7f1d1d": "Al vessel", "#4a90d9": "³He gas bore"}[col])}/>')
o5.append(sd.line(sx(-40), cy_prof, sx(56), cy_prof, "#555", 1, "6 5"))
o5.append(sd.T(sx(-36), cy_prof - 148, "nose (faces the beam)", 21, sd.MUT, "start"))
o5.append(sd.T(sx(54), cy_prof - 70, "neck / valve", 21, sd.MUT, "end"))
o5.append(sd.line(sx(-20), 258, sx(-20), 600, "#999", 1.2, "4 4")); o5.append(sd.line(sx(20), 258, sx(20), 600, "#999", 1.2, "4 4"))
o5.append(sd.T(sx(0), 270, "barrel", 19, sd.MUT))
for z in range(-30, 51, 10):
    o5.append(sd.T(sx(z), 650, f"{z}", 18, sd.MUT))
o5.append(sd.T(sx(10), 676, "z [mm]", 20, sd.MUT))
cards = sd.col(
    sd.card(sd.p(f'<b style="color:{sd.RED}">{db["nose_butt_frac"] * 100:.0f}%</b> of accepted tracks leave through the dome or neck', 26) +
            sd.p(f'{db["x_tot_nosebutt"] * 100:.2f}% X₀ of capsule there vs {db["x_tot_barrel"] * 100:.2f}% through the barrel '
                 f'(mean Al path {db["al_mm_nosebutt"]:.1f} vs {db["al_mm_barrel"]:.2f} mm).', 22, sd.MUT),
            tip="Ray tracing of 3×10⁵ isotropic tracks from vertices uniform in the gas, kept if they hit a MM arm; exit point = first step outside the CFRP."),
    sd.card(sd.p("The simulated baseline already contains them", 26) +
            sd.p("pairs_v2 used the full Geant4 geometry, and the upgrade ratio uses the same ray-traced budgets, so nothing is ignored.", 22, sd.MUT)),
    sd.card(sd.p(f'Upgrade cell: <b style="color:{sd.BLUE}">{du["frac_lost_endcap"] * 100:.1f}%</b> exit through an end cap', 26) +
            sd.p(f"Dropped from the acceptance, not counted as scattering. Skin path factor ⟨1/cos⟩ = {du['path_factor_cell']:.2f}.", 22, sd.MUT)),
    w=600, gap=16)
D.slide("nose",
        sd.title(f"A quarter of baseline tracks leave through the dome or neck, where the capsule is {db['x_tot_nosebutt'] / db['x_tot_barrel']:.1f}× thicker",
                 "Ray-traced exit position of tracks that reach a MM arm; vertices uniform in the ³He gas (pairs_v2 STEP capsule)") +
        sd.row(sd.svg(PW, PH + 20, "".join(o5), "capsule profile with exit-position histogram"), cards, gap=44),
        short="Nose / butt",
        notes="Acceptance: four flat arms at 250 mm, active area 399 mm tangential × 360 mm along the beam, centred on the origin; pinwheel shifts ignored. "
              "The 'barrel' is |z_exit| ≤ 20 mm. This answers whether the budget is a side-exit-only estimate: it is not, "
              "and the nose and neck (solid Al, r = 3.5 mm) matter because roughly a quarter of the accepted tracks cross them.")

# ── slide 6: the X17 peak ────────────────────────────────────────────────────
P6 = sd.Plot(1060, 600, x=(60, 180), y=(0, .2), xlabel="e⁺e⁻ opening angle [deg]",
             ylabel="fraction of X17 pairs per degree")
P6.xticks([(k, str(k)) for k in range(60, 181, 20)])
P6.yticks([(k / 100, f"{k / 100:.2f}") for k in range(0, 21, 5)])
xs_ = [c for c in cen if 60 <= c <= 180]
sel = [i for i, c in enumerate(cen) if 60 <= c <= 180]
P6.band(xs_, [0] * len(xs_), [ht[i] for i in sel], sd.GREEN, alpha=.18)
for key, colr, nm in (("ideal", sd.GREEN, "ideal"), ("base", sd.RED, "Baseline"), ("up", sd.BLUE, "Upgrade")):
    yv = [ht[i] if key == "ideal" else RES[key]["h"][i] for i in sel]
    P6.line(xs_, yv, colr, w=3.5, markers=False,
            tips=[f"{nm}\nθ = {a:.0f}°: {b:.4f} per degree" for a, b in zip(xs_, yv)])
P6.vline(109.4, sd.GREY, label="θ_min = 109°", tip="Kinematic shoulder of the at-rest X17 (m = 16.8 MeV): θ_min = 2 asin(m/E*), E* = 20.58 MeV.")
rows6 = [[sd.term("Baseline", "STEP capsule, standard MM window"), f'{RES["base"]["s68"]:.1f}°', f'×{RES["base"]["rel"]:.2f}', f'{RES["base"]["fwhm"]:.0f}°'],
         [sd.term("Upgrade", "2 bar mylar cell + mylar window"), f'{RES["up"]["s68"]:.1f}°', f'×{RES["up"]["rel"]:.2f}', f'{RES["up"]["fwhm"]:.0f}°'],
         ["Ideal", "0°", "×1.00", "–"]]
D.slide("peak",
        sd.title(f"The 109° shoulder returns: peak height ×{RES['base']['rel']:.2f} → ×{RES['up']['rel']:.2f}, FWHM {RES['base']['fwhm']:.0f}° → {RES['up']['fwhm']:.0f}°",
                 "Simulated X17 opening angle after smearing each leg by the rescaled Geant4 direction error; 2×10⁶ at-rest pairs, acceptance not applied") +
        sd.row(P6.svg("X17 opening-angle peak"),
               sd.col(sd.legend([("ideal", sd.GREEN), ("baseline", sd.RED), ("upgrade", sd.BLUE)], size=24),
                      sd.table(["", "σ68(Δθ)", "peak", "FWHM"], rows6, size=24, widths=[170, 110, 100, 100]),
                      sd.callout(f"Air and drift gas are {(mu['air'] + mu['drift']) / xu * 100:.0f}% of what remains.", sd.BLUE, 24),
                      sd.callout("A sharp shoulder is back. Whether cut-and-count works also depends on the IPC background, which is not modelled here.", sd.GREEN, 24),
                      w=560, gap=22), gap=44),
        short="X17 peak",
        foot="Direction at the first MM hit (no vertex constraint, valid for the larger cell); drift gas added to both cases; peak height = 3° running mean relative to the truth.",
        notes=("Method: sample X17 pairs at rest (E* = 20.58 MeV), tilt each leg by ψ drawn from P(ψ|KE) of the 'first' estimator in pairs_v2, "
               "multiplied by θ₀(x_case)/θ₀(x_baseline upstream) at that KE. This keeps the measured non-Gaussian tails. The baseline reproduces the note's "
               f"σ68 = 15.5°. Sensitivity: with a 7 µm rolled-Al permeation layer on the cell skin (instead of 2×40 nm), x/X₀ rises "
               f"{xu * 100:.3f}% → {F.xtot(m7) * 100:.3f}% and σ68 {RES['up']['s68']:.1f}° → {RES['up7']['s68']:.1f}°. "
               "Not modelled: acceptance changes (the 300 mm cell has a longer vertex distribution than the 60 mm capsule), energy loss and range of "
               "soft legs through the old wall, rods and seam of the cell (≈ 9% of azimuth), IPC background."))

# ── slide 7: appendix, 4He bag ───────────────────────────────────────────────
o7 = []
yb7, sc7 = 640, 640 / (xu * 1.0)
for xc, m, nm in ((150, mu, "Upgrade (air gap)"), (440, ma, "Upgrade + ⁴He bag")):
    s, sp, _ = stack_svg_parts(m, xc, 150, yb7, sc7, labels=False)
    o7.append(s)
    o7.append(sd.T(xc, yb7 + 34, nm, 22, sd.INK, "middle", 600))
    o7.append(sd.T(xc, yb7 - F.xtot(m) * sc7 - 12, f"{F.xtot(m) * 100:.3f}%", 24, sd.INK, "middle", 600))
for t in np.arange(0, xu * 100 + 1e-9, 0.04):
    y = yb7 - t / 100 * sc7
    o7.append(sd.line(56, y, 64, y, "#888", 1.2)); o7.append(sd.T(50, y + 7, f"{t:.2f}%", 17, sd.MUT, "end"))
# legend of components in the right-hand bar
s_, sp_, seg_ = stack_svg_parts(ma, 440, 150, yb7, sc7, labels=True, lab_x=540, lab_gap=30, size=19)
o7.append(s_)
P7 = sd.Plot(760, 560, x=(80, 180), y=(0, .2), xlabel="e⁺e⁻ opening angle [deg]", ylabel="fraction of X17 pairs per degree")
P7.xticks([(k, str(k)) for k in range(80, 181, 20)]); P7.yticks([(k / 100, f"{k / 100:.2f}") for k in range(0, 21, 5)])
xs7 = [c for c in cen if 80 <= c <= 180]; sel7 = [i for i, c in enumerate(cen) if 80 <= c <= 180]
P7.band(xs7, [0] * len(xs7), [ht[i] for i in sel7], sd.GREEN, alpha=.18)
for key, colr, nm in (("ideal", sd.GREEN, "ideal"), ("up", sd.BLUE, "Upgrade"), ("app", sd.ORANGE, "Upgrade + ⁴He bag")):
    yv = [ht[i] if key == "ideal" else RES[key]["h"][i] for i in sel7]
    P7.line(xs7, yv, colr, w=3.5, markers=False, tips=[f"{nm}\nθ = {a:.0f}°: {b:.4f} per degree" for a, b in zip(xs7, yv)])
D.slide("appendix",
        sd.title(f"Appendix: ⁴He instead of air, plus two foils, halves it again — σ68 {RES['up']['s68']:.1f}° → {RES['app']['s68']:.1f}°",
                 f"The air gap (0.084%) becomes {ma['he4'] * 100:.3f}% of ⁴He plus {(ma['bag_Mylar'] + ma['bag_Al']) * 100:.3f}% of foil; the 30 mm drift gas is then the largest term") +
        sd.row(sd.svg(1000, 720, "".join(o7), "budget of the upgrade with and without a 4He bag"),
               sd.col(P7.svg("peak with 4He bag"),
                      sd.legend([("ideal", sd.GREEN), ("upgrade", sd.BLUE), ("+ ⁴He bag", sd.ORANGE)], size=22), w=800), gap=24),
        short="Appendix: ⁴He bag",
        notes=(f"Scenario: the 24 cm between target and MM is a ⁴He volume (ρ = 0.166 mg/cm³, X₀ = 5.7×10⁵ cm), closed by two further 12 µm aluminised-mylar foils "
               f"(one at each end of the bag). Peak height ×{RES['up']['rel']:.2f} → ×{RES['app']['rel']:.2f}, FWHM {RES['up']['fwhm']:.0f}° → {RES['app']['fwhm']:.0f}°. "
               f"The remaining budget is {F.xtot(ma) * 100:.3f}% X₀, of which {ma['drift'] / F.xtot(ma) * 100:.0f}% is the drift gas, so removing the bag volume altogether "
               "(operating the MM in the target vacuum or gas volume) would be the only way to go further. The bag foils must withstand the pressure difference "
               "to air; the inner window foil holds the HV."))

# ── slide 8: limits and decisions ────────────────────────────────────────────
lim = sd.row(
    sd.card(sd.p("What this does not rule out", 30, sd.DINK, 600) +
            sd.p("• The upgrade acceptance is not modelled: a 300 mm cell has a much longer vertex distribution than the 60 mm capsule, "
                 "and low-energy legs the old wall stopped would now reach the MM.", 24, sd.DMUT) +
            sd.p("• Foil thicknesses are assumptions (12 µm window foils; 40 nm Al coatings). A 7 µm Al permeation layer on the skin costs "
                 f"σ68 {RES['up']['s68']:.1f}° → {RES['up7']['s68']:.1f}°.", 24, sd.DMUT) +
            sd.p("• Rods and seam (≈ 9% of azimuth) scatter far more than the facets and are ignored.", 24, sd.DMUT) +
            sd.p("• Scaling a measured table by the mean material ratio ignores that individual tracks cross very different thicknesses.", 24, sd.DMUT),
            bg="#1f2533", w=800, pad=36) +
    sd.card(sd.p("Decisions this points at", 30, sd.DINK, 600) +
            sd.p(f"1. Go to the low-pressure mylar cell: it removes {share_wall * 100:.0f}% of the baseline budget.", 24, sd.DMUT) +
            sd.p("2. Choose the MM window by what holds HV, not by the foil stack: two mylar foils cost 0.01% X₀.", 24, sd.DMUT) +
            sd.p(f"3. Then the air is {share_air * 100:.0f}% of what remains: decide whether a ⁴He bag (σ68 → {RES['app']['s68']:.1f}°) is worth the extra foils and engineering.", 24, sd.DMUT) +
            sd.p("4. Re-run with a full Geant4 simulation of the cell geometry to replace the scaling.", 24, sd.DMUT),
            bg="#1f2533", w=800, pad=36), gap=48)
D.slide("limits",
        sd.kicker("Caveats and next steps") +
        '<h2 style="font-size:60px;font-weight:600;line-height:1.1;letter-spacing:-1px">'
        'The scattering we had at n_TOF was large; at low pressure it becomes a question about the air</h2>' + lim,
        dark=True, short="Limits",
        notes="Provenance: angular_resolution_note.md (the baseline study this builds on), vessel_design in x17_facility_search (cell geometry), "
              "scripts/scattering_geometry.py and make_scattering_upgrade_figs.py in MX17_Full_Geant.")

D.write(OUT, note_meta=dict(
    title=f"Reducing MX17 Material Budget → {RES['base']['s68']:.0f}° to {RES['up']['s68']:.0f}° scattering",
    summary=(f"Ray-traced material budgets and X17 opening-angle peak: {xb * 100:.1f}% → {xu * 100:.2f}% X₀, "
             f"σ68 {RES['base']['s68']:.0f}° → {RES['up']['s68']:.1f}°; appendix with a ⁴He bag."),
    tags="X17, nTOF, Geant4, multiple scattering, opening angle, target",
    date="2026-10-06"))
print("wrote", OUT)
