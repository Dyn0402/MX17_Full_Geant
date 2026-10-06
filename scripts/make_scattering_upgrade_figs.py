#!/usr/bin/env python3
"""make_scattering_upgrade_figs.py — what do two material reductions buy us?

BASELINE : STEP capsule (0.6 mm Al + 0.9 mm CFRP, 500 bar) + MM window/cathode
           (40 um Mylar + 50 um Kapton + 9 um Cu)                [pairs_v2 sim]
UPGRADE  : (1) 2 bar 3He mylar-wrap cell (x17_facility_search/vessel_design)
           (2) MM drift window -> two aluminised mylar foils (gas + HV)
APPENDIX : UPGRADE with the 24 cm air gap replaced by a 4He bag closed by two
           more aluminised-mylar foils.

Material budgets are acceptance-averaged by ray tracing (scattering_geometry.py):
real capsule profile incl. nose/butt exits, wall obliquity, oblique path through
air / window / drift gas.  The 30 mm drift gas is included in every case.

Resolution: the baseline per-leg direction error P(psi | KE) is the Geant4
response ("first" estimator = true direction at the first MM hit, pairs_v2);
each psi is rescaled by the Highland ratio theta0(x_case)/theta0(x_baseline
upstream) at the same KE, keeping the measured non-Gaussian tails.  Rough
estimate: no acceptance / energy-loss / range changes, no vertex constraint.

Usage:  python scripts/make_scattering_upgrade_figs.py
"""
import sys
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Rectangle

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_angular_resolution_figs as ar    # noqa: E402
import scattering_geometry as geo            # noqa: E402

OUT = REPO / "docs/scattering_upgrade"
OUT.mkdir(parents=True, exist_ok=True)
ME = ar.ME
plt.rcParams.update({"axes.axisbelow": True, "axes.grid": False,
                     "font.size": 10})

# ── components: key -> (family, label, colour) ───────────────────────────────
# Same hue = same physical part; light = mylar, mid = CFRP/Kapton/air, dark = metal
COMP = {
    "he3":        ("³He gas", "³He gas", "#4a90d9"),
    "Al":         ("Capsule / cell wall", "Al capsule (0.6 mm barrel, dome)", "#7f1d1d"),
    "CFRP":       ("Capsule / cell wall", "CFRP wrap (0.9 mm)", "#d62728"),
    "Mylar_cell": ("Capsule / cell wall", "Mylar skin (12 µm)", "#f2a1a1"),
    "Al_cell":    ("Capsule / cell wall", "Al coating (2×40 nm)", "#7f1d1d"),
    "air":        ("Gap to MM", "Air (≈24 cm)", "#2ca02c"),
    "he4":        ("Gap to MM", "⁴He bag gas (≈24 cm)", "#8fd18f"),
    "bag_Mylar":  ("Gap to MM", "2 extra mylar foils (12 µm)", "#cdeccd"),
    "bag_Al":     ("Gap to MM", "their Al coating", "#1e6b1e"),
    "win_Cu":     ("MM window", "Cu cathode (9 µm)", "#4b2c85"),
    "win_Kapton": ("MM window", "Kapton (50 µm)", "#9467bd"),
    "win_Mylar":  ("MM window", "Mylar (40 µm)", "#d3c2ea"),
    "winf_Mylar": ("MM window", "2 window foils: mylar (12 µm each)", "#d3c2ea"),
    "winf_Al":    ("MM window", "their Al coating", "#4b2c85"),
    "drift":      ("Drift gas", "Drift gas (30 mm ArIso)", "#8c8c8c"),
}
ORDER = ["he3", "Al", "CFRP", "Al_cell", "Mylar_cell", "air", "he4", "bag_Al",
         "bag_Mylar", "win_Cu", "win_Kapton", "winf_Al", "win_Mylar", "winf_Mylar", "drift"]


def stack(m):
    return [(k, m[k]) for k in ORDER if k in m]


def xtot(m):
    return sum(m.values())


def highland(ke, x):
    p = np.sqrt((ke + ME) ** 2 - ME ** 2)
    return ar.highland(p, x)


# ── measured baseline median psi(KE) ("first" estimator) ─────────────────────
def med_sim(d):
    de = d["direction_error"]
    ke_e = np.array(de["ke_bin_edges"]); ke_c = .5 * (ke_e[:-1] + ke_e[1:])
    psi_e = np.array(de["psi_bin_edges"]); psi_c = .5 * (psi_e[:-1] + psi_e[1:])
    h = np.array(de["estimators"]["first"]["counts"], float)
    med = np.array([ar.hist_quantile(r, psi_c, .5) for r in h])
    ok = ~np.isnan(med)
    return lambda ke: np.interp(ke, ke_c[ok], med[ok])


def build_cases():
    xb, db = geo.baseline()
    xu, du = geo.upgrade()
    xa, da = geo.upgrade(he4_gap=True, n_extra_foils=2)
    ren = lambda m: {k.replace("win_", "winf_"): v for k, v in m.items()}
    mb, mu, ma = geo.mean(xb), ren(geo.mean(xu)), ren(geo.mean(xa))
    return dict(base=mb, up=mu, app=ma), dict(base=db, up=du, app=da)


# ── Fig 1 ────────────────────────────────────────────────────────────────────
def draw_stack(ax, x0, w, m, unit=100, label_min=None, fmt="{:.3f}%"):
    bot = 0.0
    for k, v in stack(m):
        fam, lab, col = COMP[k]
        ax.bar(x0, v * unit, bottom=bot, width=w, color=col, edgecolor="white",
               lw=.6)
        bot += v * unit
    return bot


def legend_by_family(ax, keys, **kw):
    from matplotlib.patches import Patch
    handles, labels, heads, last = [], [], [], None
    for k in [k for k in ORDER if k in keys][::-1]:
        fam, lab, col = COMP[k]
        if fam != last:
            handles.append(Patch(color="none")); labels.append(fam); heads.append(len(labels) - 1)
            last = fam
        handles.append(Patch(color=col)); labels.append("  " + lab)
    leg = ax.legend(handles, labels, frameon=False, fontsize=8, handlelength=1.2,
                    handletextpad=.5, labelspacing=.25, **kw)
    for i in heads:
        leg.get_texts()[i].set_fontweight("bold")
    return leg


def fig1(M, med):
    mb, mu = M["base"], M["up"]
    xb, xu = xtot(mb), xtot(mu)
    xb_up = xb - mb["drift"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 6.4),
                                 gridspec_kw=dict(width_ratios=[1.25, 1]))
    tb = draw_stack(a1, 0, .55, mb)
    tu = draw_stack(a1, 1, .55, mu)
    a1.text(0, tb + .03, f"{tb:.2f}%", ha="center", fontweight="bold")
    a1.text(1, tu + .03, f"{tu:.2f}%", ha="center", fontweight="bold")
    a1.set_xticks([0, 1])
    a1.set_xticklabels(["Baseline\n(500 bar capsule,\nKapton/Cu window)",
                        "Upgrade\n(2 bar mylar cell,\nmylar window)"])
    a1.set_ylabel("material before the MM strips, x / X₀  [%]")
    a1.set_title(f"Material the e⁺e⁻ cross (acceptance-averaged)  —  "
                 f"{xb / xu:.0f}× less", fontsize=11)
    a1.set_xlim(-.45, 1.45); a1.set_ylim(0, 3.8)
    a1.grid(axis="y", alpha=.25)

    # exploded view of the upgrade stack, in the empty space above it
    ins = a1.inset_axes([.52, .22, .45, .70])
    zmax = tu * 1.08
    segs, bot = [], 0.0
    for k, v in stack(mu):
        fam, lab, col = COMP[k]
        ins.bar(0, v * 100, bottom=bot, width=.5, color=col, edgecolor="white", lw=.6)
        segs.append([bot + v * 50, lab, v * 100, k]); bot += v * 100
    # spread the labels so they do not collide, then draw leader lines
    ys = [sg[0] for sg in segs]; gap = zmax * .075
    for _ in range(200):
        for i in range(1, len(ys)):
            if ys[i] - ys[i - 1] < gap:
                m = .5 * (ys[i] + ys[i - 1]); ys[i - 1] = m - gap / 2; ys[i] = m + gap / 2
    if ys[0] < gap * .6:                       # keep labels inside the frame
        ys = [y + (gap * .6 - ys[0]) for y in ys]
        for _ in range(200):
            for i in range(1, len(ys)):
                if ys[i] - ys[i - 1] < gap:
                    ys[i] = ys[i - 1] + gap
    for (yb, lab, v, k), yl in zip(segs, ys):
        ins.plot([.25, .45], [yb, yl], color="#777", lw=.6)
        ins.text(.47, yl, f"{lab}   {v:.4f}%" if v < .01 else f"{lab}   {v:.3f}%",
                 va="center", fontsize=7.5)
    ins.set_xlim(-.4, 3.5); ins.set_ylim(0, zmax)
    ins.set_xticks([]); ins.set_ylabel("x / X₀ [%]  (zoomed)", fontsize=8)
    ins.tick_params(labelsize=8)
    ins.set_title("Upgrade stack, exploded", fontsize=9)
    ins.grid(axis="y", alpha=.25)
    for sp in ins.spines.values():
        sp.set_color("#555")
    for xy_a, xy_b in (((.725, tu), (0, 0)), ((1.275, tu), (1, 0))):
        a1.add_artist(ConnectionPatch(xyA=xy_a, coordsA=a1.transData, xyB=xy_b,
                                      coordsB=ins.transAxes, color="#777", lw=.8, ls="--"))
    a1.add_patch(Rectangle((.725, 0), .55, tu, fill=False, ec="#777", lw=.8))
    legend_by_family(a1, set(mb) | set(mu), loc="upper left",
                     bbox_to_anchor=(0.0, 1.0))

    # right: median scattering angle from the measured baseline, scaled
    kes = [3, 5, 10]
    w = .34
    for j, (m, col, lab) in enumerate([(mb, "#d62728", "Baseline"),
                                       (mu, "#1f77b4", "Upgrade")]):
        vals = [med(k) * highland(k, xtot(m)) / highland(k, xb_up) for k in kes]
        bars = a2.bar(np.arange(3) + (j - .5) * w, vals, w, color=col, label=lab)
        for b_, v in zip(bars, vals):
            a2.text(b_.get_x() + w / 2, v + .3, f"{v:.1f}°", ha="center", fontsize=9)
    a2.set_xticks(range(3)); a2.set_xticklabels([f"{k} MeV" for k in kes])
    a2.set_xlabel("e± kinetic energy"); a2.set_ylabel("median scattering angle per track [deg]")
    a2.set_title("Per-track scattering at the MM entrance\n"
                 "(X17 legs: KE⁻+KE⁺ ≈ 19.6 MeV → soft leg ≤ 9.8 MeV)", fontsize=11)
    a2.legend(frameon=False); a2.grid(axis="y", alpha=.25)
    fig.tight_layout(); fig.savefig(OUT / "fig1_budget_bars.png", dpi=160)
    plt.close(fig)


# ── Fig 2 ────────────────────────────────────────────────────────────────────
def fig2(M, med):
    mb, mu = M["base"], M["up"]
    xb_up = xtot(mb) - mb["drift"]
    ke = np.linspace(2, 15, 200)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), sharey=True)
    for ax, m, ttl in zip(axes, (mb, mu), ("Baseline", "Upgrade")):
        xt = xtot(m)
        tot = med(ke) * highland(ke, xt) / highland(ke, xb_up)
        bot = np.zeros_like(ke)
        # stack in order, band height = share of the scattering variance
        for k, v in stack(m):
            fam, lab, col = COMP[k]
            ax.fill_between(ke, bot, bot + v / xt * tot, color=col, lw=0,
                            label=f"{lab} ({v / xt * 100:.0f}%)" if v / xt > .005 else f"{lab} (<1%)")
            bot += v / xt * tot
        ax.plot(ke, tot, "k", lw=1.4)
        ax.set_title(f"{ttl}: x/X₀ = {xt * 100:.2f}%")
        ax.set_xlabel("e± kinetic energy [MeV]"); ax.grid(alpha=.25)
        h, l = ax.get_legend_handles_labels()
        ax.legend(h[::-1], l[::-1], frameon=False, fontsize=8, loc="upper right")
    axes[0].set_ylabel("median scattering angle per track [deg]")
    axes[0].set_ylim(0, 36); axes[0].set_xlim(2, 15)
    fig.suptitle("Where the scattering comes from vs e± energy — band height = "
                 "share of the scattering variance (∝ x/X₀)")
    fig.text(.5, .005, "Baseline total = measured Geant4 median (direction at first MM "
             "hit); upgrade = same, rescaled by the Highland ratio.", ha="center",
             fontsize=8, color="#555")
    fig.tight_layout(rect=(0, .02, 1, 1)); fig.savefig(OUT / "fig2_scattering_vs_energy.png", dpi=160)
    plt.close(fig)


# ── Fig 3 / A1: X17 peak ─────────────────────────────────────────────────────
def sigma68(d):
    lo, hi = np.percentile(d, [16, 84])
    return .5 * (hi - lo)


def peak_curves(M, d, cases, n=2_000_000, seed=7):
    rng = np.random.default_rng(seed)
    th, ke1, ke2, u1, u2 = ar.sample_pairs(ar.M_X17, n, rng)
    base_s = ar.make_psi_sampler(d, "first")
    xb_up = xtot(M["base"]) - M["base"]["drift"]
    bins = np.linspace(0, 180, 181)
    cen = .5 * (bins[:-1] + bins[1:])
    sm = lambda h, w=3: np.convolve(h, np.ones(w) / w, mode="same")
    ht = np.histogram(th, bins=bins)[0] / n
    res = {}
    for key in cases:
        x = xtot(M[key])
        ratio = lambda ke, x=x: highland(ke, x) / highland(ke, xb_up)
        samp = lambda ke, r, ratio=ratio: base_s(ke, r) * ratio(ke)
        thr = ar.smear_pair_response(ke1, ke2, u1, u2, samp, rng)
        h = np.histogram(thr, bins=bins)[0] / n
        hs = sm(h)
        half = cen[hs >= hs.max() / 2]
        res[key] = dict(h=hs, s68=sigma68(thr - th), rel=hs.max() / sm(ht).max(),
                        fwhm=half.max() - half.min())
    return cen, ht, res


def draw_peak(ax, cen, ht, res, spec):
    ax.fill_between(cen, ht, step="mid", color="#2ca02c", alpha=.2)
    ax.step(cen, ht, where="mid", color="#2ca02c", lw=2, label="ideal (no scattering)")
    for key, (lbl, col, ls) in spec.items():
        r = res[key]
        ax.step(cen, r["h"], where="mid", color=col, ls=ls, lw=2,
                label=f"{lbl}:  σ68 = {r['s68']:.1f}°, peak ×{r['rel']:.2f}")
    ax.set_xlim(60, 180); ax.set_ylim(0, .2)
    ax.set_xlabel("e⁺e⁻ opening angle [deg]"); ax.set_ylabel("fraction of X17 pairs per degree")
    ax.grid(alpha=.25); ax.legend(frameon=False, fontsize=9, loc="upper right")


def peak_bars(ax, res, spec):
    labs = [v[0].replace(" ", "\n", 1) for v in spec.values()]
    cols = [v[1] for v in spec.values()]
    vals = [res[k]["rel"] for k in spec]
    bars = ax.bar(labs, vals, color=cols)
    ax.axhline(1, color="#2ca02c", lw=2)
    ax.text(len(vals) - .55, 1.02, "ideal", color="#2ca02c", ha="right", fontsize=9)
    for b_, k in zip(bars, spec):
        ax.text(b_.get_x() + b_.get_width() / 2, res[k]["rel"] + .02,
                f"{res[k]['rel']:.2f}\nFWHM {res[k]['fwhm']:.0f}°", ha="center", fontsize=9)
    ax.set_ylim(0, 1.25); ax.set_ylabel("peak height relative to ideal")
    ax.set_title("Peak sharpness"); ax.grid(axis="y", alpha=.25)


def fig3(M, d):
    spec = {"base": ("Baseline", "#d62728", "-"), "up": ("Upgrade", "#1f77b4", "-")}
    cen, ht, res = peak_curves(M, d, spec)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw=dict(width_ratios=[1.6, 1]))
    draw_peak(axes[0], cen, ht, res, spec)
    axes[0].set_title("Expected X17 opening-angle peak (at-rest capture)\n"
                      "direction at first MM hit, drift gas included; acceptance not applied")
    peak_bars(axes[1], res, spec)
    fig.tight_layout(); fig.savefig(OUT / "fig3_x17_peak.png", dpi=160)
    plt.close(fig)
    return res


def figA1(M, d):
    spec = {"up": ("Upgrade (air gap)", "#1f77b4", "-"),
            "app": ("Upgrade + ⁴He bag", "#ff7f0e", "-")}
    cen, ht, res = peak_curves(M, d, spec)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.4),
                             gridspec_kw=dict(width_ratios=[1.0, 1.5, .8]))
    ax = axes[0]
    for i, k in enumerate(("up", "app")):
        draw_stack(ax, i, .55, M[k])
        ax.text(i, xtot(M[k]) * 100 + .005, f"{xtot(M[k]) * 100:.3f}%", ha="center", fontweight="bold")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Upgrade\n(air gap)", "Upgrade\n+ ⁴He bag"])
    ax.set_ylabel("x / X₀ [%]"); ax.set_title("Material budget"); ax.grid(axis="y", alpha=.25)
    ax.set_xlim(-.5, 2.9); ax.set_ylim(0, .17)
    legend_by_family(ax, set(M["up"]) | set(M["app"]), loc="upper right", bbox_to_anchor=(1.02, 1.0))
    draw_peak(axes[1], cen, ht, res, spec)
    axes[1].set_title("X17 peak")
    peak_bars(axes[2], res, spec)
    fig.tight_layout(); fig.savefig(OUT / "figA1_he4_bag.png", dpi=160)
    plt.close(fig)
    return res


def main():
    d = json.load(open(ar.RESP))
    med = med_sim(d)
    M, diag = build_cases()
    for k, m in M.items():
        print(k, f"x/X0 = {xtot(m) * 100:.3f}%", {a: round(b * 100, 4) for a, b in m.items()})
    print(json.dumps(diag, indent=1))
    fig1(M, med); fig2(M, med)
    r3 = fig3(M, d); rA = figA1(M, d)
    for nm, r in (("main", r3), ("appendix", rA)):
        for k, v in r.items():
            print(f"{nm:9s}{k:5s} sigma68={v['s68']:.2f} peak=x{v['rel']:.2f} FWHM={v['fwhm']:.0f}")
    json.dump(dict(M=M, diag=diag, main={k: {a: b for a, b in v.items() if a != 'h'} for k, v in r3.items()},
                   appendix={k: {a: b for a, b in v.items() if a != 'h'} for k, v in rA.items()}),
              open(OUT / "results.json", "w"), indent=1)


if __name__ == "__main__":
    main()
