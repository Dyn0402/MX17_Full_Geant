#!/usr/bin/env python3
"""
make_thermal_report_figures.py
Generate all figures for docs/report/thermal_note.tex from a
analyze_thermal_captures.py .npz output.

Usage (lxplus):
    source scripts/setup_lxplus.sh
    python3 scripts/make_thermal_report_figures.py thermal_captures_subkev_prelim.npz \
        --outdir docs/report/figs [--lambda2d data/lamda2DvsEn_EAR2.root]
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

plt.rcParams.update({"font.size": 11, "figure.dpi": 120,
                     "axes.grid": True, "grid.alpha": 0.3})

SIGMA_RATIO = 54e-6 / 5333.0          # sigma_ng / sigma_np at thermal
ALBERTO = {1e-2: 1.41e-6, 1e-1: 5.71e-7, 1e0: 1.68e-7,
           1e1: 6.69e-8, 1e2: 5.90e-8}
CAP = {1e-2: 1.01e-8, 1e-1: 1.01e-8, 1e0: 1.01e-8, 1e1: 1.00e-8, 1e2: 7.9e-9}
DEPTH_SLICES = [(1e-3, 1e-1), (1e-1, 1e1), (1e1, 1e3)]


def poisson_interval(n):
    """Approximate 68% central Poisson interval (Garwood-style)."""
    if n == 0:
        return 0.0, 1.84
    lo = n * (1 - 1/(9*n) - 1/(3*np.sqrt(n)))**3
    hi = (n+1) * (1 - 1/(9*(n+1)) + 1/(3*np.sqrt(n+1)))**3
    return lo, hi


def decade_sums(h):
    return np.array([h[i*10:(i+1)*10].sum() for i in range(6)])


# ── fig 1: reaction tree ──────────────────────────────────────────────────────
def fig_reaction(out):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)

    def box(x, y, w, h, text, fc="#eef3fb", fontsize=10.5):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12",
                                    fc=fc, ec="k", lw=1.1))
        ax.text(x + w/2, y + h/2, text, ha="center", va="center",
                fontsize=fontsize)

    def arrow(x1, y1, x2, y2, label="", dx=0.15, fs=9.5, color="k"):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=16, lw=1.4, color=color))
        if label:
            ax.text((x1+x2)/2 + dx, (y1+y2)/2, label, fontsize=fs, color=color)

    box(0.4, 7.6, 3.6, 1.6, "slow neutron + $^3$He\n(the gas in the capsule)")
    box(5.8, 7.6, 3.6, 1.6, "$^4$He$^*$\nexcited by 20.6 MeV")
    arrow(4.0, 8.4, 5.8, 8.4, "absorbed", dx=-0.4)

    box(0.4, 4.4, 3.8, 1.5, "proton + triton\n(heat in the gas, no photon)",
        fc="#fbeee6")
    arrow(6.6, 7.6, 3.0, 5.9)
    ax.text(3.6, 6.9, "almost always\n($\\sigma$ = 5333 b at 25 meV)",
            fontsize=9.5, ha="right")

    box(5.8, 4.4, 3.8, 1.5,
        "$^4$He + 20.6 MeV de-excitation\n(the only path to X17)",
        fc="#e8f6e8")
    arrow(7.7, 7.6, 7.7, 5.9,
          "1 in $10^8$\n($\\sigma$ = 54 $\\mu$b)", dx=0.15)

    box(3.6, 0.6, 2.4, 1.4, "real $\\gamma$\n(most of the time)")
    box(0.4, 0.6, 2.4, 1.4, "e$^+$e$^-$ pair (IPC)\n2.1 per 1000")
    box(6.8, 0.6, 2.8, 1.4, "X17 $\\to$ e$^+$e$^-$\n(if it exists)",
        fc="#fff3cc")
    arrow(7.0, 4.4, 4.9, 2.0)
    arrow(6.9, 4.4, 1.9, 2.0)
    arrow(8.2, 4.4, 8.2, 2.0)

    ax.set_title("What can happen when a slow neutron meets $^3$He",
                 fontsize=13)
    fig.savefig(out / "fig_reaction_tree.pdf", bbox_inches="tight")
    plt.close(fig)


# ── fig 2: transparent vs opaque cartoon (real capsule, beam from below) ─────
def _draw_capsule(ax):
    """Side view of the STEP capsule [mm]: gas bore r=10, 40 mm barrel +
    hemispherical ends (60 mm on axis); Al dome below, neck+valve above."""
    ax.set_xlim(-45, 45); ax.set_ylim(-80, 66)
    ax.set_aspect("equal"); ax.axis("off")
    # Aluminium silhouette (slightly schematic dome/neck)
    t = np.linspace(np.pi, 2 * np.pi, 40)
    dome = np.column_stack([10.6 * np.cos(t), -20 + 15 * np.sin(t)])
    al = np.vstack([[[-10.6, 20]], dome[::-1], [[10.6, 20]],
                    [[3.5, 27]], [[3.5, 51]], [[-3.5, 51]], [[-3.5, 27]]])
    ax.add_patch(plt.Polygon(al, closed=True, fc="#9a9a9a", ec="k",
                             lw=0.8, zorder=2))
    # Gas: 40 mm cylinder + 10 mm hemispherical ends
    t1 = np.linspace(np.pi, 2 * np.pi, 40)
    t2 = np.linspace(0, np.pi, 40)
    gas = np.vstack([
        np.column_stack([10 * np.cos(t1), -20 + 10 * np.sin(t1)]),
        np.column_stack([10 * np.cos(t2), 20 + 10 * np.sin(t2)]),
    ])
    ax.add_patch(plt.Polygon(gas, closed=True, fc="#dcebf7", ec="k",
                             lw=0.8, zorder=3))
    # beam arrows from below
    for x in (-5.0, 0.0, 5.0):
        ax.add_patch(FancyArrowPatch((x, -58), (x, -38), arrowstyle="-|>",
                                     mutation_scale=13, color="C0", lw=1.6,
                                     zorder=4))
    ax.text(13, -50, "beam", color="C0", fontsize=9)
    # 60 mm dimension marker
    ax.annotate("", xy=(20, 30), xytext=(20, -30),
                arrowprops=dict(arrowstyle="<->", color="0.4", lw=1.1))
    ax.text(22, 0, "60 mm of gas\non the beam axis", fontsize=8.5,
            color="0.3", va="center")
    ax.text(-13, -33, "5 mm\nAl dome", fontsize=8, ha="right", color="0.3")


def fig_opaque(out):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 6.0))

    ax = axes[0]
    _draw_capsule(ax)
    ax.set_title("The assumption: gas is transparent\n"
                 "(thin-target formula)", fontsize=11)
    for x in (-5.0, 0.0, 5.0):  # neutrons sail through the whole column
        ax.add_patch(FancyArrowPatch((x, -28), (x, 34), arrowstyle="-|>",
                                     mutation_scale=12, color="C0", lw=1.3,
                                     ls=(0, (4, 2)), zorder=5))
    ax.text(0, -62, "every neutron samples the whole 60 mm column:\n"
            "P(rare) = (atoms/cm$^2$) $\\times\\ \\sigma_{n\\gamma}$\n"
            "grows with thickness", ha="center", va="top", fontsize=9)

    ax = axes[1]
    _draw_capsule(ax)
    ax.set_title("The reality below 1 keV: gas is opaque\n"
                 "(mean free path $\\approx$ 0.15 mm at 25 meV)", fontsize=11)
    # absorption band just inside the gas entrance (thickness exaggerated)
    tb = np.linspace(np.pi + 0.25, 2 * np.pi - 0.25, 40)
    band = np.vstack([
        np.column_stack([10.0 * np.cos(tb), -20 + 10.0 * np.sin(tb)]),
        np.column_stack([7.5 * np.cos(tb[::-1]), -20 + 7.5 * np.sin(tb[::-1])]),
    ])
    ax.add_patch(plt.Polygon(band, closed=True, fc="#c23b22", alpha=0.8,
                             zorder=4))
    for x in (-5.0, 0.0, 5.0):
        ax.plot([x], [-20 - 10 * np.sqrt(1 - (x / 10) ** 2) + 1.5], marker="*",
                ms=12, color="#7a1f12", zorder=6)
    ax.annotate("neutrons never reach\nthe rest of the gas",
                xy=(-3, 0), xytext=(-43, 8), fontsize=9, ha="left",
                va="center",
                arrowprops=dict(arrowstyle="->", color="0.3", lw=1.1,
                                shrinkB=4))
    ax.text(0, -62, "absorbed within a fraction of a mm of entering\n"
            "(red band exaggerated);  P(rare) $\\to\\ \\sigma_{n\\gamma}/"
            "\\sigma_{np} = 1.0\\times10^{-8}$\nthickness-independent",
            ha="center", va="top", fontsize=9)

    fig.savefig(out / "fig_opaque_cartoon.pdf", bbox_inches="tight")
    plt.close(fig)


# ── fig 3: beam per pulse + where neutrons end up vs E ───────────────────────
def fig_beam_fates(d, out):
    edges = d["loge_edges"]; centers = 10.0 ** ((edges[:-1] + edges[1:]) / 2)
    w = float(d["n_per_pulse"]) / float(d["n_events"])
    h_all = d["h_all"].astype(float)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    ax1.step(centers, h_all * w, where="mid", color="k", lw=1.5)
    ax1.set_yscale("log"); ax1.set_ylabel("neutrons / pulse / bin")
    ax1.set_title("EAR2 beam below 1 keV as sampled by the simulation\n"
                  "(evaluated flux $\\times$ energy-dependent footprint)")

    with np.errstate(divide="ignore", invalid="ignore"):
        for h, lab, c in [(d["h_gas_np"],  "absorbed in the gas (n,p)t", "C0"),
                          (d["h_wall_al"], "captured in Al capsule",     "C1"),
                          (d["h_scint_h"], "captured in scintillators",  "C2")]:
            frac = np.where(h_all > 0, h.astype(float) / h_all, np.nan)
            ax2.step(centers, frac, where="mid", label=lab, color=c, lw=1.5)
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel("neutron energy [eV]")
    ax2.set_ylabel("fraction of beam neutrons")
    ax2.legend(fontsize=9); ax2.set_ylim(1e-4, 2)
    fig.savefig(out / "fig_beam_fates.pdf", bbox_inches="tight")
    plt.close(fig)


# ── fig 4: money plot ─────────────────────────────────────────────────────────
def fig_money(d, out):
    edges = d["loge_edges"]
    w = float(d["n_per_pulse"]) / float(d["n_events"])
    np_dec   = decade_sums(d["h_gas_np"].astype(float))
    rad_dec  = decade_sums(d["h_gas_rad"].astype(float))
    beam_dec = decade_sums(d["h_all"].astype(float))
    dec_edges = 10.0 ** np.arange(-3, 4)
    dec_cent  = np.sqrt(dec_edges[:-1] * dec_edges[1:])

    fig, ax = plt.subplots(figsize=(9, 6.3))
    ax.step(dec_edges, np.append(np_dec, np_dec[-1]) * SIGMA_RATIO * w,
            where="post", color="C0", lw=2.2,
            label="full simulation: (n,p)t counts $\\times\\ "
                  "\\sigma_{n\\gamma}/\\sigma_{np}$")
    for i in range(6):
        n = rad_dec[i]
        if n > 0:
            lo, hi = poisson_interval(n)
            ax.errorbar([dec_cent[i]], [n * w],
                        yerr=[[(n - lo) * w], [(hi - n) * w]],
                        fmt="o", color="C3", ms=9, capsize=5, lw=1.8, zorder=6)
        elif beam_dec[i] > 0:
            ul = 1.84 * w
            ax.annotate("", xy=(dec_cent[i], ul * 0.4),
                        xytext=(dec_cent[i], ul),
                        arrowprops=dict(arrowstyle="-|>", color="C3", lw=1.6))
    ax.errorbar([], [], fmt="o", color="C3", ms=9, capsize=5,
                label="full simulation: direct $^3$He(n,$\\gamma$) counts")

    tt_x, tt_y, cp_y, sim_y = [], [], [], []
    for i, elo in enumerate(10.0 ** np.arange(-3.0, 3.0)):
        for key in ALBERTO:
            if abs(np.log10(key / elo)) < 0.05:
                tt_x.append(dec_cent[i])
                tt_y.append(ALBERTO[key] * beam_dec[i] * w)
                cp_y.append(CAP[key] * beam_dec[i] * w)
                sim_y.append(np_dec[i] * SIGMA_RATIO * w)
    ax.plot(tt_x, tt_y, "s--", color="C1", ms=8,
            label="thin-target rate table (transparent-gas assumption)")
    ax.plot(tt_x, cp_y, "^:", color="C2", ms=8,
            label="opaque-gas ceiling $\\sigma_{n\\gamma}/\\sigma_{np}$ "
                  "$\\times$ beam")
    # table / simulation overshoot factors (match Appendix A)
    for x, ty, sy in zip(tt_x, tt_y, sim_y):
        ax.annotate(f"$\\times${ty/sy:.0f}", xy=(x, np.sqrt(ty * sy)),
                    fontsize=9, color="0.35", ha="left",
                    xytext=(x * 1.25, np.sqrt(ty * sy)))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("neutron energy [eV]")
    ax.set_ylabel("$^3$He(n,$\\gamma$)$^4$He per pulse per energy decade")
    tot = np_dec.sum() * SIGMA_RATIO * w
    ax.set_title("X17-production captures per pulse, sub-keV region\n"
                 f"simulation total: {tot:.1e} /pulse  "
                 f"($\\to$ {tot*2.1e-3:.1e} IPC/pulse;  "
                 "table: 1.2$\\times$10$^{-2}$ IPC/pulse)")
    ax.legend(fontsize=9, loc="upper right")
    fig.savefig(out / "fig_money.pdf", bbox_inches="tight")
    plt.close(fig)


# ── fig 5: absorption depth ───────────────────────────────────────────────────
def fig_depth(d, out):
    capy_edges = d["capy_edges"]
    capy_c = (capy_edges[:-1] + capy_edges[1:]) / 2
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = [f"{lo:g}–{hi:g} eV" for (lo, hi) in DEPTH_SLICES]
    for i, lab in enumerate(labels):
        h = d["h_depth"][i].astype(float)
        if h.sum() > 0:
            ax.step(capy_c, h / h.sum() / (capy_edges[1] - capy_edges[0]),
                    where="mid", lw=1.7, label=lab)
    # mark the 4 direct radiative events
    rad = d["rad_events"]
    if rad.size:
        for E, y in rad:
            ax.axvline(y, color="#c23b22", alpha=0.6, lw=1.2, ls="--")
        ax.plot([], [], color="#c23b22", ls="--",
                label="the 4 direct (n,$\\gamma$) events")
    ax.set_xlabel("where along the beam axis the neutron was absorbed [mm]\n"
                  "(beam enters from the left)")
    ax.set_ylabel("probability density [mm$^{-1}$]")
    ax.set_yscale("log")
    ax.set_title("Self-shielding seen directly: absorption depth vs "
                 "neutron energy")
    ax.legend(fontsize=9)
    fig.savefig(out / "fig_depth.pdf", bbox_inches="tight")
    plt.close(fig)


# ── fig 6: footprint vs bore (from Lambda2D if available) ────────────────────
def fig_footprint(out, lambda2d_path):
    fig, ax = plt.subplots(figsize=(8, 5))
    drew_data = False
    if lambda2d_path and Path(lambda2d_path).is_file():
        try:
            import uproot
            with uproot.open(lambda2d_path) as f:
                h2 = f["Lambda2D"]
                vals = h2.values()                       # (x, y)
                xe = h2.axis(0).edges(); ye = h2.axis(1).edges()
                # Identify the energy axis: the one spanning many decades
                if xe[0] > 0 and np.log10(xe[-1] / max(xe[0], 1e-30)) > 4:
                    e_edges, r_edges, m = xe, ye, vals
                else:
                    e_edges, r_edges, m = ye, xe, vals.T
                e_cent = np.sqrt(e_edges[:-1] * np.maximum(e_edges[1:], 1e-30))
                r_cent = (r_edges[:-1] + r_edges[1:]) / 2
                r90 = np.full(len(e_cent), np.nan)
                for i in range(len(e_cent)):
                    row = m[i]
                    if row.sum() > 0:
                        c = np.cumsum(row) / row.sum()
                        r90[i] = np.interp(0.9, c, r_cent)
                ok = ~np.isnan(r90)
                ax.plot(e_cent[ok], r90[ok], lw=1.8, color="C0",
                        label="radius containing 90% of the beam (EAR2 data)")
                ax.set_xscale("log")
                drew_data = True
        except Exception as e:
            print(f"  (Lambda2D read failed: {e}; drawing schematic only)")
    ax.axhline(1.0, color="#c23b22", lw=2,
               label="capsule bore radius (10 mm)")
    ax.set_xlabel("neutron energy [eV]")
    ax.set_ylabel("beam radius [cm]")
    ax.set_title("Beam footprint vs capsule bore: the halo misses the gas")
    ax.legend(fontsize=9)
    if not drew_data:
        ax.text(0.5, 0.55, "(EAR2 Lambda2D file not available\n"
                "when this figure was generated)",
                transform=ax.transAxes, ha="center")
    fig.savefig(out / "fig_footprint.pdf", bbox_inches="tight")
    plt.close(fig)


# ── fig 5b: measured absorption positions inside the capsule ─────────────────
# The "cartoon, but real": 2D map of where beam neutrons are absorbed in the
# gas, thermal slice vs keV slice, drawn inside the true capsule profile.
GAS_Z  = np.array([-30., -28., -26., -24., -22., -20.,
                    20.,  22.,  24.,  26.,  28.,  30.])          # mm
GAS_RO = np.array([1e-3, 6.0, 8.0, 9.165, 9.798, 10.0,
                   10.0, 9.798, 9.165, 8.0, 6.0, 1e-3])
AL_Z   = np.array([-35.000, -34.635, -33.947, -32.880, -31.447,
                   -29.693, -27.688, -25.521, -23.288, -21.075,
                   -20.000,  20.000,  22.000,  24.000,  26.000,
                    26.770,  31.370,  40.900,  50.980])
AL_RO  = np.array([0.000, 2.323, 3.902, 5.433, 6.850,
                   8.090, 9.102, 9.856, 10.342, 10.573,
                   10.600, 10.600, 10.317, 9.469, 8.054,
                   7.360, 6.912, 3.500, 3.500])


def _capsule_outline(ax):
    for z, ro, lw, c in [(AL_Z, AL_RO, 1.3, "0.25"),
                         (GAS_Z, GAS_RO, 0.9, "0.45")]:
        ax.plot(np.concatenate([ro, -ro[::-1], [ro[0]]]),
                np.concatenate([z, z[::-1], [z[0]]]), color=c, lw=lw,
                zorder=6)


def fig_gas_absorption(out, root_file):
    import uproot
    slices = [("thermal:  $E_n$ < 0.1 eV", 0.0, 0.1),
              ("0.1–1 keV", 100.0, 1000.0)]
    xs = [[] for _ in slices]; ys = [[] for _ in slices]
    with uproot.open(root_file) as f:
        for chunk in f["EventTree"].iterate(
                ["neutron_E_eV", "capture_vol", "capture_proc",
                 "cap_x", "cap_y"],
                step_size=2_000_000, library="np"):
            vol  = np.asarray(chunk["capture_vol"],  dtype=object)
            proc = np.asarray(chunk["capture_proc"], dtype=object)
            gas = (vol == "He3Gas") & (proc == "neutronInelastic")
            E = chunk["neutron_E_eV"]
            for i, (_, elo, ehi) in enumerate(slices):
                m = gas & (E >= elo) & (E < ehi)
                xs[i].append(chunk["cap_x"][m])
                ys[i].append(chunk["cap_y"][m])

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 6.2), sharey=True)
    for i, (ax, (lab, _, _)) in enumerate(zip(axes, slices)):
        x = np.concatenate(xs[i]); y = np.concatenate(ys[i])
        hb = ax.hexbin(x, y, gridsize=110,
                       extent=(-12, 12, -36, 36), bins="log",
                       cmap="inferno", mincnt=1)
        _capsule_outline(ax)
        ax.add_patch(FancyArrowPatch((0, -52), (0, -38), arrowstyle="-|>",
                                     mutation_scale=14, color="C0", lw=1.8))
        ax.text(2.5, -47, "beam", color="C0", fontsize=10)
        ax.set_xlim(-14, 14); ax.set_ylim(-55, 38)
        ax.set_aspect("equal")
        ax.set_xlabel("x [mm]")
        ax.set_title(f"{lab}\n({len(x):,} absorptions)", fontsize=11)
        ax.grid(False)
        fig.colorbar(hb, ax=ax, shrink=0.8, label="absorptions (log)")
    axes[0].set_ylabel("y along beam [mm]")
    fig.suptitle("Where beam neutrons are actually absorbed in the gas "
                 "(simulation, one job file)", fontsize=12.5)
    fig.savefig(out / "fig_gas_absorption.pdf", bbox_inches="tight")
    plt.close(fig)


# ── fig 6b: where scintillator-captured neutrons end up (from one run-B file) ─
def fig_scint_origin(out, root_file):
    import uproot
    SCINT = ("PlasticScint", "LiqScint_1", "LiqScint_2",
             "BackScintL", "BackScintR")
    xs, ys, zs = [], [], []
    with uproot.open(root_file) as f:
        for chunk in f["EventTree"].iterate(
                ["capture_vol", "capture_proc", "cap_x", "cap_y", "cap_z"],
                step_size=2_000_000, library="np"):
            vol  = np.asarray(chunk["capture_vol"],  dtype=object)
            proc = np.asarray(chunk["capture_proc"], dtype=object)
            m = np.isin(vol, SCINT) & (proc == "nCapture")
            xs.append(chunk["cap_x"][m]); ys.append(chunk["cap_y"][m])
            zs.append(chunk["cap_z"][m])
    x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5),
                                   gridspec_kw={"width_ratios": [1.15, 1]})
    hb = ax1.hexbin(x / 10, z / 10, gridsize=80, bins="log", cmap="viridis")
    fig.colorbar(hb, ax=ax1, label="captures (log scale)")
    ax1.add_patch(plt.Circle((0, 0), 1.15, fc="none", ec="r", lw=1.5))
    ax1.annotate("capsule", xy=(1.2, 1.2), color="r", fontsize=9)
    ax1.set_aspect("equal")
    ax1.set_xlabel("x [cm]"); ax1.set_ylabel("z [cm]")
    ax1.set_title("Capture positions, beam view\n(four-arm pattern)")

    ax2.hist(y / 10, bins=80, histtype="step", lw=1.6)
    ax2.set_xlabel("y along beam [cm]"); ax2.set_ylabel("captures / bin")
    ax2.set_title("Capture positions along the beam axis")
    fig.suptitle(f"Neutron captures in the scintillator volumes "
                 f"({len(x):,} in one job file)", fontsize=12)
    fig.savefig(out / "fig_scint_origin.pdf", bbox_inches="tight")
    plt.close(fig)


# ── fig 7: normalization ladder ───────────────────────────────────────────────
def fig_ladder(d, out):
    w = float(d["n_per_pulse"]) / float(d["n_events"])
    n_np = float(d["h_gas_np"].sum())
    rad_pp = n_np * SIGMA_RATIO * w
    steps = [
        ("EAR2 beam < 1 keV", "$7.31\\times10^{6}$ n/pulse", "#eef3fb"),
        ("absorbed in the $^3$He gas",
         f"$\\times$ {n_np/float(d['n_events']):.2f}  $\\to$  "
         f"{n_np*w:.2e} /pulse", "#eef3fb"),
        ("radiative branch (n,$\\gamma$)",
         f"$\\times\\ 1.0\\times10^{{-8}}$  $\\to$  {rad_pp:.1e} /pulse",
         "#e8f6e8"),
        ("internal e$^+$e$^-$ pairs (IPC)",
         f"$\\times\\ 2.1\\times10^{{-3}}$  $\\to$  "
         f"{rad_pp*2.1e-3:.1e} /pulse", "#fbeee6"),
        ("X17 (assumed 2.5%)",
         f"$\\times\\ 0.025$  $\\to$  {rad_pp*0.025:.1e} /pulse", "#fff3cc"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 6.6))
    ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(-0.5, 10.4)
    y = 9.0
    for title, val, fc in steps:
        ax.add_patch(FancyBboxPatch((1.2, y - 0.65), 7.6, 1.3,
                                    boxstyle="round,pad=0.1", fc=fc, ec="k"))
        ax.text(5.0, y + 0.22, title, ha="center", fontsize=11)
        ax.text(5.0, y - 0.30, val, ha="center", fontsize=10.5)
        if y > 1.5:
            ax.add_patch(FancyArrowPatch((5.0, y - 0.78), (5.0, y - 1.32),
                                         arrowstyle="-|>", mutation_scale=18,
                                         lw=1.5))
        y -= 2.1
    ax.set_title("Per-pulse normalisation ladder (sub-keV, full simulation)",
                 fontsize=12.5)
    fig.savefig(out / "fig_ladder.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--outdir", default="docs/report/figs")
    ap.add_argument("--lambda2d", default="data/lamda2DvsEn_EAR2.root")
    ap.add_argument("--scint-events", default=None, metavar="ROOT_FILE",
                    help="One run-B ROOT file; makes the scintillator-capture"
                         " position figure (slow: full EventTree scan)")
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    d = np.load(args.npz, allow_pickle=True)

    fig_reaction(out)
    fig_opaque(out)
    fig_beam_fates(d, out)
    fig_money(d, out)
    fig_depth(d, out)
    fig_footprint(out, args.lambda2d)
    fig_ladder(d, out)
    if args.scint_events:
        fig_gas_absorption(out, args.scint_events)
        fig_scint_origin(out, args.scint_events)
    print(f"Figures written to {out}/")


if __name__ == "__main__":
    main()
