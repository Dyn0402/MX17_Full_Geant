#!/usr/bin/env python3
"""make_angular_report_figures.py — figures for docs/report/angular_note.tex

Produces (into docs/report/figs/):
  fig_ang_diagram.pdf   summary scattering diagram: top-down real geometry,
                        one worked X17 event, zoom inset on the capsule wall
  fig_ang_gallery.pdf   3-panel gallery: same truth pair, three scatter outcomes
  fig_ang_fan.pdf       scatter cones (68/95%) on the geometry + the resulting
                        opening-angle distribution from a 3D toy MC
  fig_psi_vs_ke.pdf     |
  fig_sigma68_vs_minke.pdf | PDF versions of the analysis figures from
  fig_ms_budget.pdf     |   make_angular_resolution_figs.py (same code)
  fig_theta_dilution.pdf|
  fig_ipc_dilution.pdf  |  IPC truth vs smeared spectrum
  fig_theta_money.pdf   |  money plot: X17 truth vs smeared, alone + on IPC

Geometry is the current SimConfig/DetectorConstruction (STEP capsule,
MM front face at 250 mm — note plot_geometry.py still carries 22 cm).

Usage:  python scripts/make_angular_report_figures.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Arc, Wedge, FancyArrowPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_angular_resolution_figs as mar

REPO = Path(__file__).resolve().parent.parent
OUT  = REPO / "docs/report/figs"

# ── Geometry (cm), top-down XZ plane at y = 0  (SimConfig current) ───────────
DIST      = 25.0          # MM window front face from origin
R_GAS     = 1.00          # He-3 bore radius
R_AL      = 1.06          # Al barrel outer radius
R_CFRP    = 1.15          # CFRP wrap outer radius
R_SCAT    = 1.10          # effective scatter radius (mid-wall) used in cartoons

# arm stack offsets from the MM front face [cm] (thicknesses from
# DetectorConstruction.cc; thin foils lumped into the drift-gas front edge)
W_DRIFT   = (0.01, 3.01)          # 30 mm drift gas
W_PCB     = (3.05, 3.61)
W_TRIG    = (5.61, 5.95)          # after 20 mm air gap
W_LS      = (7.95, 12.68)         # LS stack incl. CFRP walls
W_BSC     = (12.78, 15.22)        # back plastic bars
HW_MM, HW_TRIG, HW_LS, HW_BSC = 19.0, 24.0, 22.5, 25.4

ARM_DIRS = {0: np.array([1, 0]), 1: np.array([-1, 0]),
            2: np.array([0, 1]), 3: np.array([0, -1])}

# ── Worked example kinematics ────────────────────────────────────────────────
THETA_TRUE = 109.4        # symmetric at-rest X17 decay opening angle [deg]
KE_SYM     = (mar.E_X17 - 2 * mar.ME) / 2          # 9.78 MeV per leg
P_SYM      = np.sqrt((KE_SYM + mar.ME)**2 - mar.ME**2)
XX0_TOT    = sum(t / x0 for _, t, x0 in mar.MS_LAYERS)
THETA0_SYM = np.degrees(mar.highland(P_SYM, XX0_TOT))   # ≈ 7.1° plane-projected

BIS = 45.0                                          # decay bisector angle [deg]
ANG_EM = BIS - THETA_TRUE / 2                       # e- toward +X arm
ANG_EP = BIS + THETA_TRUE / 2                       # e+ toward +Z arm


def _u(deg):
    r = np.radians(deg)
    return np.array([np.cos(r), np.sin(r)])


def draw_arms(ax, lw=0.5, alpha=0.8, full_stack=True):
    """Draw the four detector arms (top-down XZ)."""
    layers = [(W_DRIFT, HW_MM, "#4a90d9", "MM drift gas (30 mm)"),
              (W_PCB,   HW_MM, "#5cb85c", "MM PCB")]
    if full_stack:
        layers += [(W_TRIG, HW_TRIG, "#f0c040", "trigger plastic scint."),
                   (W_LS,   HW_LS,   "#d9534f", "liquid-scint. stack"),
                   (W_BSC,  HW_BSC,  "#e07820", "back plastic scints")]
    seen = set()
    for aid, w_hat in ARM_DIRS.items():
        u_hat = np.array([-w_hat[1], w_hat[0]])
        for (w0, w1), hu, col, lab in layers:
            cen = w_hat * (DIST + (w0 + w1) / 2)
            hw, hus = (w1 - w0) / 2, hu
            hx = abs(w_hat[0]) * hw + abs(u_hat[0]) * hus
            hz = abs(w_hat[1]) * hw + abs(u_hat[1]) * hus
            ax.add_patch(Rectangle((cen[0] - hx, cen[1] - hz), 2 * hx, 2 * hz,
                                   lw=lw, edgecolor="k", facecolor=col,
                                   alpha=alpha, zorder=2,
                                   label=lab if lab not in seen else None))
            seen.add(lab)


def draw_capsule(ax, lw=0.6):
    for r, fc, zo, lab in [(R_CFRP, "#404040", 3, "CFRP wrap (0.9 mm)"),
                           (R_AL,   "#b0b0b0", 4, "Al barrel (0.6 mm)"),
                           (R_GAS,  "#99d8f5", 5, "He-3 gas bore (r = 10 mm)")]:
        ax.add_patch(plt.Circle((0, 0), r, fc=fc, ec="k", lw=lw,
                                zorder=zo, label=lab))


def kinked_track(vtx, ang_true, psi, t_end):
    """Vertex → wall scatter point → straight continuation.
    Returns (scatter_point, end_point, ang_out).  psi = in-plane kink [deg]
    (positive rotates counter-clockwise)."""
    d0 = _u(ang_true)
    # distance from vtx to radius R_SCAT along d0
    b = vtx @ d0
    s = -b + np.sqrt(b * b + R_SCAT**2 - vtx @ vtx)
    p_sc = vtx + s * d0
    ang_out = ang_true + psi
    p_end = p_sc + t_end * _u(ang_out)
    return p_sc, p_end, ang_out


def _angle_arc(ax, cen, r, a1, a2, color, lw=1.6, ls="-", zorder=8):
    ax.add_patch(Arc(cen, 2 * r, 2 * r, angle=0, theta1=min(a1, a2),
                     theta2=max(a1, a2), color=color, lw=lw, ls=ls,
                     zorder=zorder))


def _event(ax, vtx, psi_em, psi_ep, t_em=27.5, t_ep=27.5,
           show_unscattered=True, show_backproj=True, lw=2.0):
    """Draw one pair event; returns (theta_reco, scatter points, out angles)."""
    res = {}
    for tag, ang, psi, t_end, col in [("em", ANG_EM, psi_em, t_em, "#c81e1e"),
                                      ("ep", ANG_EP, psi_ep, t_ep, "#1e50c8")]:
        p_sc, p_end, a_out = kinked_track(vtx, ang, psi, t_end)
        if show_unscattered:
            far = vtx + (np.linalg.norm(p_end - vtx) + 1.0) * _u(ang)
            ax.plot([vtx[0], far[0]], [vtx[1], far[1]], color="#2ca02c",
                    lw=1.2, ls="--", zorder=6)
        ax.plot([vtx[0], p_sc[0]], [vtx[1], p_sc[1]], color=col, lw=lw, zorder=7)
        ax.plot([p_sc[0], p_end[0]], [p_sc[1], p_end[1]], color=col, lw=lw,
                zorder=7)
        if show_backproj:
            back = p_sc - 6.0 * _u(a_out)
            ax.plot([p_sc[0], back[0]], [p_sc[1], back[1]], color=col,
                    lw=1.1, ls=":", zorder=6)
        res[tag] = (p_sc, p_end, a_out)
    ax.plot(*vtx, marker="*", ms=12, color="k", zorder=9)
    theta_reco = res["ep"][2] - res["em"][2]
    return theta_reco, res


# ─────────────────────────────────────────────────────────────────────────────
def fig_diagram():
    fig, ax = plt.subplots(figsize=(10.5, 10.5))
    ax.set_aspect("equal")
    draw_capsule(ax)
    draw_arms(ax)

    vtx = np.array([0.0, 0.0])
    psi_em, psi_ep = -9.0, +7.0          # opens the pair: Δθ = +16°
    th_reco, res = _event(ax, vtx, psi_em, psi_ep, lw=2.2)

    # truth / reco angle arcs at the vertex
    _angle_arc(ax, (0, 0), 4.2, ANG_EM, ANG_EP, "#2ca02c")
    ax.annotate(f"$\\theta_{{true}} = {THETA_TRUE:.0f}^\\circ$",
                xy=4.2 * _u(15), xytext=(12, 6.5), fontsize=12.5,
                color="#1a7a1a", zorder=9,
                arrowprops=dict(arrowstyle="->", color="#1a7a1a", lw=1.2),
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
    _angle_arc(ax, (0, 0), 6.6, res["em"][2], res["ep"][2], "#444444", ls="--")
    ax.annotate(f"$\\theta_{{reco}} = {th_reco:.0f}^\\circ$ (MM directions)",
                xy=6.6 * _u(65), xytext=(-0.5, 15.5), fontsize=12.5,
                color="#222222", zorder=9,
                arrowprops=dict(arrowstyle="->", color="#222222", lw=1.2),
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))

    # displacement at the MM between unscattered and actual hit
    for tag, arm_w in [("em", np.array([1, 0])), ("ep", np.array([0, 1]))]:
        p_sc, p_end, a_out = res[tag]
        ang = ANG_EM if tag == "em" else ANG_EP
        t = (DIST - p_sc @ arm_w) / (_u(a_out) @ arm_w)
        hit = p_sc + t * _u(a_out)
        t0 = DIST / (_u(ang) @ arm_w)
        hit0 = vtx + t0 * _u(ang)
        ax.annotate("", xy=hit, xytext=hit0,
                    arrowprops=dict(arrowstyle="<->", color="#666666", lw=1.3))
        mid = 0.5 * (hit + hit0) + np.array([0.6, 0.6])
        ax.text(*mid, f"{np.linalg.norm(hit - hit0):.0f} cm", fontsize=10,
                color="#444444")

    # labels
    ax.annotate("scatter in the capsule wall\n(~11 mm from the vertex)",
                xy=res["em"][0], xytext=(9, -13), fontsize=11,
                arrowprops=dict(arrowstyle="->", color="k", lw=1.1),
                ha="center", zorder=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
    ax.annotate("first measured point:\nMM drift gas, 250 mm away",
                xy=res["em"][1] - 1.5 * _u(res["em"][2]), xytext=(15, -22),
                fontsize=11, ha="center", zorder=9,
                arrowprops=dict(arrowstyle="->", color="k", lw=1.1),
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
    ax.plot([], [], color="#2ca02c", ls="--", lw=1.4,
            label="true directions (no scatter)")
    ax.plot([], [], color="#c81e1e", lw=2.2, label="e$^-$ actual path")
    ax.plot([], [], color="#1e50c8", lw=2.2, label="e$^+$ actual path")
    ax.plot([], [], color="k", ls=":", lw=1.2,
            label="measured track, back-projected")

    # ── zoom inset on the e- wall crossing ───────────────────────────────
    axi = ax.inset_axes([0.015, 0.015, 0.385, 0.385])
    axi.set_aspect("equal")
    draw_capsule(axi, lw=1.0)
    p_sc = res["em"][0]
    d_in, d_out = _u(ANG_EM), _u(res["em"][2])
    axi.plot([vtx[0], p_sc[0]], [vtx[1], p_sc[1]], color="#c81e1e", lw=2.4,
             zorder=7)
    axi.plot([p_sc[0], (p_sc + 1.6 * d_out)[0]],
             [p_sc[1], (p_sc + 1.6 * d_out)[1]], color="#c81e1e", lw=2.4,
             zorder=7)
    axi.plot([vtx[0], (vtx + 2.6 * d_in)[0]], [vtx[1], (vtx + 2.6 * d_in)[1]],
             color="#2ca02c", lw=1.4, ls="--", zorder=6)
    back = p_sc - 2.4 * d_out
    axi.plot([(p_sc + 1.6 * d_out)[0], back[0]],
             [(p_sc + 1.6 * d_out)[1], back[1]], color="k", ls=":", lw=1.4,
             zorder=8)
    axi.plot(*vtx, marker="*", ms=14, color="k", zorder=9)
    _angle_arc(axi, tuple(p_sc), 1.05, ANG_EM, res["em"][2], "#c81e1e", lw=1.8)
    axi.annotate(f"$\\psi = {abs(psi_em):.0f}^\\circ$ kink\nin 1.5 mm of wall",
                 xy=p_sc + 1.0 * _u(0.5 * (ANG_EM + res["em"][2])),
                 xytext=(1.7, -2.45), fontsize=10, ha="center",
                 arrowprops=dict(arrowstyle="->", lw=1.0))
    # back-projection miss distance at the vertex
    miss = np.linalg.norm(np.cross(np.append(d_out, 0),
                                   np.append(vtx - p_sc, 0)))[()]
    axi.annotate("back-projected track misses the\n"
                 f"vertex by only {10*miss:.1f} mm —\n"
                 "still points at the target!",
                 xy=back + 0.4 * d_out, xytext=(-1.15, 2.05), fontsize=9.5,
                 ha="center", arrowprops=dict(arrowstyle="->", lw=1.0))
    axi.set_xlim(-3.0, 3.0); axi.set_ylim(-3.0, 3.0)
    axi.set_xticks([]); axi.set_yticks([])
    for s in axi.spines.values():
        s.set_linewidth(1.6)
    ax.indicate_inset_zoom(axi, edgecolor="k", lw=1.4)

    ax.legend(loc="upper right", fontsize=10, framealpha=0.92)
    lim = DIST + 16.5
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("X  [cm]"); ax.set_ylabel("Z  [cm]")
    ax.set_title("How the capsule wall destroys the opening angle — "
                 "top-down view (beam into page)\n"
                 f"symmetric X17 pair, KE = {KE_SYM:.1f} MeV per leg; "
                 f"wall kinks of {abs(psi_em):.0f}$^\\circ$/{psi_ep:.0f}"
                 f"$^\\circ$ ($\\approx1\\sigma$) turn "
                 f"{THETA_TRUE:.0f}$^\\circ$ into {th_reco:.0f}$^\\circ$")
    ax.grid(True, lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(OUT / "fig_ang_diagram.pdf")
    fig.savefig(OUT / "fig_ang_diagram.png", dpi=140)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
def fig_gallery():
    cases = [(-7.0, -5.0, "closing scatter ($\\approx1\\sigma$ inward)"),
             (+7.0, +5.0, "opening scatter ($\\approx1\\sigma$ outward)"),
             (+14.0, +11.0, "tail scatter ($\\approx2\\sigma$)")]
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.0))
    for ax, (p_em, p_ep, lab) in zip(axes, cases):
        ax.set_aspect("equal")
        draw_capsule(ax)
        draw_arms(ax, full_stack=False)
        th_reco, res = _event(ax, np.array([0.0, 0.0]), -p_em, p_ep,
                              show_backproj=False, lw=1.8)
        _angle_arc(ax, (0, 0), 4.0, ANG_EM, ANG_EP, "#2ca02c")
        _angle_arc(ax, (0, 0), 6.0, res["em"][2], res["ep"][2], "#444444",
                   ls="--")
        ax.set_title(f"{lab}\n$\\theta_{{true}} = {THETA_TRUE:.0f}^\\circ"
                     f"\\;\\rightarrow\\;\\theta_{{reco}} = {th_reco:.0f}"
                     f"^\\circ$  ($\\Delta\\theta = "
                     f"{th_reco-THETA_TRUE:+.0f}^\\circ$)", fontsize=11)
        lim = DIST + 5
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_xlabel("X  [cm]")
        ax.grid(True, lw=0.3, alpha=0.4)
    axes[0].set_ylabel("Z  [cm]")
    fig.suptitle("The same true pair, three throws of the scattering dice "
                 f"(per-leg $\\theta_0 \\approx {THETA0_SYM:.0f}^\\circ$ at "
                 f"KE = {KE_SYM:.1f} MeV)", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "fig_ang_gallery.pdf")
    fig.savefig(OUT / "fig_ang_gallery.png", dpi=140)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
def fig_fan(n=400_000, seed=7):
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.6),
                             gridspec_kw={"width_ratios": [1.05, 1]})

    # left: scatter cones on the geometry
    ax = axes[0]
    ax.set_aspect("equal")
    draw_capsule(ax)
    draw_arms(ax, full_stack=False)
    vtx = np.array([0.0, 0.0])
    for ang, col in [(ANG_EM, "#c81e1e"), (ANG_EP, "#1e50c8")]:
        p_sc, _, _ = kinked_track(vtx, ang, 0.0, 1.0)
        ax.plot([vtx[0], p_sc[0]], [vtx[1], p_sc[1]], color=col, lw=2.0,
                zorder=7)
        for half, alpha in [(2 * THETA0_SYM, 0.18), (THETA0_SYM, 0.38)]:
            ax.add_patch(Wedge(tuple(p_sc), 28.5, ang - half, ang + half,
                               color=col, alpha=alpha, zorder=4))
        far = vtx + 29.0 * _u(ang)
        ax.plot([vtx[0], far[0]], [vtx[1], far[1]], color=col, lw=1.0,
                ls="--", zorder=6)
    ax.plot(*vtx, marker="*", ms=12, color="k", zorder=9)
    _angle_arc(ax, (0, 0), 4.0, ANG_EM, ANG_EP, "#2ca02c")
    ax.text(*(5.2 * _u(BIS)), f"$\\theta_{{true}} = {THETA_TRUE:.0f}^\\circ$",
            color="#2ca02c", fontsize=11, ha="left", va="center")
    ax.text(13, -12, f"$\\pm\\theta_0$ ({THETA0_SYM:.0f}$^\\circ$) and "
                    f"$\\pm2\\theta_0$ cones\nafter the capsule wall",
            fontsize=11, ha="center", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
    lim = DIST + 5
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("X  [cm]"); ax.set_ylabel("Z  [cm]")
    ax.set_title("Where each track can land after the wall\n"
                 "(symmetric pair, 68%/95% scattering cones)")
    ax.grid(True, lw=0.3, alpha=0.4)

    # right: 3D toy MC of the reconstructed opening angle
    ax = axes[1]
    rng = np.random.default_rng(seed)
    th0 = np.radians(THETA0_SYM)
    a = np.radians(THETA_TRUE)
    u1 = np.array([1.0, 0.0, 0.0])
    u2 = np.array([np.cos(a), np.sin(a), 0.0])

    def smear(u):
        e1 = np.cross(u, [0, 0, 1.0]); e1 /= np.linalg.norm(e1)
        e2 = np.cross(u, e1)
        t1 = rng.normal(0, th0, n); t2 = rng.normal(0, th0, n)
        v = u[None, :] + t1[:, None] * e1[None, :] + t2[:, None] * e2[None, :]
        return v / np.linalg.norm(v, axis=1, keepdims=True)

    v1, v2 = smear(u1), smear(u2)
    th = np.degrees(np.arccos(np.clip((v1 * v2).sum(1), -1, 1)))
    q16, q50, q84 = np.percentile(th, [16, 50, 84])
    s68 = 0.5 * (q84 - q16)
    bins = np.linspace(60, 170, 111)
    ax.hist(th, bins=bins, color="#888888", alpha=0.65, density=True,
            label=f"toy MC: two legs, $\\theta_0={THETA0_SYM:.1f}^\\circ$ each\n"
                  f"$\\sigma_{{68}} = {s68:.1f}^\\circ$, "
                  f"median bias ${q50-THETA_TRUE:+.1f}^\\circ$")
    ax.axvline(THETA_TRUE, color="#2ca02c", lw=2.2,
               label=f"$\\theta_{{true}} = {THETA_TRUE:.0f}^\\circ$")
    for q in (q16, q84):
        ax.axvline(q, color="k", lw=1.0, ls=":")
    ax.annotate("", xy=(q84, 0.028), xytext=(q16, 0.028),
                arrowprops=dict(arrowstyle="<->", color="k", lw=1.2))
    ax.text(q50, 0.0292, "68% of events", ha="center", fontsize=10)
    ax.set_xlabel("reconstructed opening angle  [deg]")
    ax.set_ylabel("probability density")
    ax.set_title("What one fixed true angle becomes\n"
                 "(Geant4 full-spectrum value: $\\sigma_{68}=12.5$–$15.5^\\circ$)")
    ax.legend(fontsize=9.5); ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUT / "fig_ang_fan.pdf")
    fig.savefig(OUT / "fig_ang_fan.png", dpi=140)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig_diagram()
    fig_gallery()
    fig_fan()
    d = json.load(open(mar.RESP))
    mar.fig_psi_vs_ke(d, out=OUT, ext="pdf")
    mar.fig_sigma68_vs_minke(d, out=OUT, ext="pdf")
    mar.fig_ms_budget(out=OUT, ext="pdf")
    mar.fig_theta_dilution(d, out=OUT, ext="pdf")
    mar.fig_ipc_dilution(out=OUT, ext="pdf")
    mar.fig_res_vs_theta(d, out=OUT, ext="pdf")
    mar.fig_theta_money(d, out=OUT, ext="pdf")
    print(f"report figures written to {OUT}")


if __name__ == "__main__":
    main()
