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
    arrow(6.6, 7.6, 3.0, 5.9,
          "almost always\n($\\sigma$ = 5333 b at 25 meV)", dx=0.25)

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


# ── fig 2: transparent vs opaque cartoon ─────────────────────────────────────
def fig_opaque(out):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
    for ax in axes:
        ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
        # gas slab
        ax.add_patch(Rectangle((3.0, 1.5), 4.0, 7.0, fc="#dcebf7", ec="k"))
        ax.text(5.0, 9.0, "30 mm of $^3$He gas", ha="center", fontsize=10)

    ax = axes[0]
    ax.set_title("The assumption: gas is transparent\n"
                 "(thin-target formula)", fontsize=11)
    for y in (3.0, 5.0, 7.0):
        ax.add_patch(FancyArrowPatch((0.5, y), (9.5, y), arrowstyle="-|>",
                                     mutation_scale=14, color="C0", lw=1.6))
    ax.text(5.0, 0.5, "every neutron samples the whole column:\n"
            "P(rare capture) = (atoms/cm$^2$) $\\times\\ \\sigma_{n\\gamma}$"
            "  — grows with thickness", ha="center", fontsize=9.5)

    ax = axes[1]
    ax.set_title("The reality below 1 keV: gas is opaque\n"
                 "(mean free path $\\approx$ 0.1 mm at 25 meV)", fontsize=11)
    ax.add_patch(Rectangle((3.0, 1.5), 0.35, 7.0, fc="#c23b22", alpha=0.75))
    for y in (3.0, 5.0, 7.0):
        ax.add_patch(FancyArrowPatch((0.5, y), (3.25, y), arrowstyle="-|>",
                                     mutation_scale=14, color="C0", lw=1.6))
        ax.plot([3.28], [y], marker="*", ms=13, color="#c23b22", zorder=5)
    ax.text(5.6, 5.0, "neutrons never\nreach the back\nof the gas",
            ha="center", fontsize=10)
    ax.text(5.0, 0.5, "every neutron is absorbed regardless; the rare branch\n"
            "can only win its share:  P(rare) $\\to\\ \\sigma_{n\\gamma}/"
            "\\sigma_{np} = 1.0\\times10^{-8}$  — thickness-independent",
            ha="center", fontsize=9.5)

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

    tt_x, tt_y, cp_y = [], [], []
    for i, elo in enumerate(10.0 ** np.arange(-3.0, 3.0)):
        for key in ALBERTO:
            if abs(np.log10(key / elo)) < 0.05:
                tt_x.append(dec_cent[i])
                tt_y.append(ALBERTO[key] * beam_dec[i] * w)
                cp_y.append(CAP[key] * beam_dec[i] * w)
    ax.plot(tt_x, tt_y, "s--", color="C1", ms=8,
            label="thin-target rate table (transparent-gas assumption)")
    ax.plot(tt_x, cp_y, "^:", color="C2", ms=8,
            label="opaque-gas ceiling $\\sigma_{n\\gamma}/\\sigma_{np}$ "
                  "$\\times$ beam")
    # factor annotations
    for x, ty, cy in zip(tt_x, tt_y, cp_y):
        ax.annotate(f"$\\times${ty/cy:.0f}", xy=(x, np.sqrt(ty * cy)),
                    fontsize=9, color="0.35", ha="left",
                    xytext=(x * 1.25, np.sqrt(ty * cy)))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("neutron energy [eV]")
    ax.set_ylabel("$^3$He(n,$\\gamma$)$^4$He per pulse per energy decade")
    tot = np_dec.sum() * SIGMA_RATIO * w
    ax.set_title("X17-production captures per pulse, sub-keV region\n"
                 f"simulation total: {tot:.1e} /pulse  "
                 f"($\\to$ {tot*2.1e-3:.1e} IPC/pulse;  "
                 "table: 1.2$\\times$10$^{-2}$ IPC/pulse)")
    ax.legend(fontsize=9, loc="upper left")
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
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10.4)
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
    print(f"Figures written to {out}/")


if __name__ == "__main__":
    main()
