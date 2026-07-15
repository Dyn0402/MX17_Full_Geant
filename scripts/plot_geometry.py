#!/usr/bin/env python3
"""
plot_geometry.py
Visualise the MX17 4-arm detector geometry.
  Figure 1 — top-down 2D (+Z right, +X up, beam +Y ⊙ out of page)  [matplotlib]
  Figure 1b — side view (XY plane)                                 [matplotlib]
  Figure 2 — 3D isometric view                                     [pyvista/VTK]

No Geant4 dependency; dimensions are taken directly from SimConfig defaults.
Run from anywhere:  python scripts/plot_geometry.py

Stack (inside → out, measured 2026-07-15; see GEOMETRY_CHANGE_CHECKLIST.md):
  MM → SiPM wall → plastics → 1 liquid-scintillator layer.
  • SiPM wall (50×50 cm) = 20 bars of 2.5 cm, centered on the STRUCTURE; only
    16 are read out (window shifted 1 bar toward the MM) — the 4 un-read bars
    are drawn transparent here and omitted from the Geant sim.
  • Plastics (2× 20×30 cm) and the single LS layer are centered on the MM
    (they inherit the per-arm pinwheel shift).
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

try:                       # pyvista only needed for the 3D view
    import pyvista as pv
    pv.OFF_SCREEN = True   # headless rendering — no display required
except ImportError:
    pv = None


# ─────────────────────────────────────────────────────────────────────────────
# SimConfig defaults  (keep in sync with include/SimConfig.hh)
# ─────────────────────────────────────────────────────────────────────────────
CFG = dict(
    # MM window front face from origin — per axis (beam at mylar-box centre):
    #   ±X arms (D,B) 40.8 cm span → 20.40;  ±Z arms (A,C) 40.9 cm span → 20.45.
    mm_distance_x_cm      = 20.40,  # ±X arms (D,B)  (SimConfig.hh)
    mm_distance_z_cm      = 20.45,  # ±Z arms (A,C)  (SimConfig.hh)
    # Per-MM tangential pinwheel shift [cm], arm order 0=D(+X) 1=B(−X) 2=A(+Z) 3=C(−Z):
    mm_pinwheel_shift_cm  = (1.55, 1.575, 1.635, 1.73),
    mm_size_u_cm          = 38.0,
    mm_size_v_cm          = 34.0,
    # SiPM trigger wall (50×50 cm; 20 bars of 2.5 cm), centered on the STRUCTURE:
    sipm_front_from_mylar_cm = 11.0,   # mylar front → SiPM container front
    sipm_container_depth_cm  =  3.3,   # container depth (scint centered inside)
    sipm_bar_width_cm        =  2.5,
    sipm_n_bars              = 20,
    sipm_n_readout           = 16,
    sipm_readout_shift_bars  =  1,     # read-out window shift toward the MM [bars]
    sipm_scint_thick_cm      =  0.3,
    sipm_size_v_cm           = 50.0,
    # Plastics (2 bars per arm, centered on the MM):
    backscint_u_cm        = 20.0,
    backscint_v_cm        = 30.0,
    backscint_thick_cm    =  2.5,   # measured ~2.5 cm (nominal 2.0)
    backscint_gap_cm      =  0.3,
    backscint_tape_um     = 200.0,
    backscint_al_um       =  20.0,
    gap_sipm_to_plastic_cm =  7.0,  # SiPM container back → plastics front (MEASURE LATER)
    # Liquid scintillator — single layer, 45×45 cm, centered on the MM:
    ls_size_u_cm          = 45.0,
    ls_size_v_cm          = 45.0,
    ls_cfrp_mm            =  2.0,
    ls_inner_cfrp_um      = 600.0,
    ls_inner_al_um        =  40.0,
    ls_thick_cm           =  2.0,
    gap_plastic_to_ls_cm  =  5.0,   # plastics back → liquid front (MEASURE LATER)
)

# ─────────────────────────────────────────────────────────────────────────────
# Layer thicknesses  [cm]
# ─────────────────────────────────────────────────────────────────────────────
_um = 1e-4
_mm = 0.1

# MM stack
tMylar    =  40 * _um
tAlWin    = 0.1 * _um
tKapCath  =  50 * _um
tCuCath   =   9 * _um
tDrift    =  30 * _mm
tMesh     =  30 * _um
tAmp      = 150 * _um
tResPaste = 100 * _um
# PCB stack
tPCBKap  =  50 * _um
tPCBCu   =  26 * _um
tPCBFR4  = 100 * _um
tPCBRoh  =   5 * _mm
tPCBAl   =  50 * _um

t_MM  = tMylar + tAlWin + tKapCath + tCuCath + tDrift + tMesh + tAmp + tResPaste
t_PCB = tPCBKap + 4*(tPCBCu + tPCBFR4) + tPCBRoh + tPCBAl

# SiPM trigger wall
bw          = CFG['sipm_bar_width_cm']
sipm_front  = CFG['sipm_front_from_mylar_cm']
sipm_contD  = CFG['sipm_container_depth_cm']
sipm_scintT = CFG['sipm_scint_thick_cm']
w_sipm_f    = sipm_front                    # container front
w_sipm_b    = sipm_front + sipm_contD       # container back
w_sipm_sc_f = sipm_front + sipm_contD/2 - sipm_scintT/2   # active scint front (thin)
w_sipm_sc_b = sipm_front + sipm_contD/2 + sipm_scintT/2   # active scint back
sipm_wall_hu = CFG['sipm_n_bars'] * bw / 2.0             # full wall half-width (u)
# Read-out window: MM sits at +t, so the window shifts toward +t (higher index).
_Nb, _Nr, _sh = CFG['sipm_n_bars'], CFG['sipm_n_readout'], CFG['sipm_readout_shift_bars']
_bc = (_Nb - 1) / 2.0 + _sh
_barLo = int(round(_bc - (_Nr - 1) / 2.0))
_barHi = int(round(_bc + (_Nr - 1) / 2.0))
SIPM_READOUT = set(range(_barLo, _barHi + 1))   # instrumented bar indices (0..Nb-1)

# Plastics (tape → Al → PVT)
tTape  = CFG['backscint_tape_um'] * _um
tBscAl = CFG['backscint_al_um']   * _um
bsc_u  = CFG['backscint_u_cm']
bsc_v  = CFG['backscint_v_cm']
bsc_th = CFG['backscint_thick_cm']
bsc_gap = CFG['backscint_gap_cm']
bscTape_hu = (bsc_u  + 2*tBscAl + 2*tTape) / 2
bscTape_hv = (bsc_v  + 2*tBscAl + 2*tTape) / 2
bscTape_hw = (bsc_th + 2*tBscAl + 2*tTape) / 2
t_bsc   = 2 * bscTape_hw
w_bsc_f = w_sipm_b + CFG['gap_sipm_to_plastic_cm']
w_bsc_b = w_bsc_f + t_bsc

# LS box (single layer): front CFRP wall | inner CFRP + Al liner | LAB | rear CFRP wall
tLSCfrp      = CFG['ls_cfrp_mm']        * _mm
tLSInnerCfrp = CFG['ls_inner_cfrp_um']  * _um
tLSInnerAl   = CFG['ls_inner_al_um']    * _um
tLS          = CFG['ls_thick_cm']       * 1.0
w_LS_f       = w_bsc_b + CFG['gap_plastic_to_ls_cm']
_lsDivFront_f = w_LS_f
_lsDivFront_b = w_LS_f + tLSCfrp + tLSInnerCfrp + tLSInnerAl   # front wall + liners
_wLS1_f       = _lsDivFront_b
_wLS1_b       = _wLS1_f + tLS
_lsDivRear_f  = _wLS1_b
_lsDivRear_b  = _wLS1_b + tLSCfrp                              # rear wall
w_LS_b        = _lsDivRear_b
t_LS_box      = w_LS_b - w_LS_f

# MM + PCB depths (from mylar front)
w_MM_f  = 0.0;        w_MM_b  = t_MM
w_PCB_f = w_MM_b;     w_PCB_b = w_PCB_f + t_PCB
stack_depth = w_LS_b

# ─────────────────────────────────────────────────────────────────────────────
# Palette
# ─────────────────────────────────────────────────────────────────────────────
C_MM, C_PCB, C_SC = '#4a90d9', '#5cb85c', '#f0c040'
C_LS, C_CFRP, C_BS = '#d9534f', '#303030', '#e07820'

# MM + PCB slabs (drawn centered on the MM)
MM_PCB_2D = [
    ('MM',  w_MM_f,  w_MM_b,  C_MM),
    ('PCB', w_PCB_f, w_PCB_b, C_PCB),
]
# LS box sub-layers (drawn centered on the MM)
LS_SUBS_2D = [
    ('LS_CFRP',   _lsDivFront_f, _lsDivFront_b, C_CFRP),
    ('LiqScint',  _wLS1_f,       _wLS1_b,       C_LS),
    ('LS_CFRP',   _lsDivRear_f,  _lsDivRear_b,  C_CFRP),
]

# ─────────────────────────────────────────────────────────────────────────────
# He-3 capsule — full STEP-derived polycone profiles [cm]
# (axis along beam = Y; keep in sync with src/DetectorConstruction.cc)
# ─────────────────────────────────────────────────────────────────────────────
Z_GAS  = np.array([-2.9500, -2.8000, -2.6000, -2.4000, -2.2000, -2.0000,
                   -1.5000, -0.5000,  0.5000,  1.5000,  2.0000,  2.2000,
                    2.4000,  2.6000,  2.8000,  3.0000,  3.2000,  3.4000,
                    3.6000,  3.8000,  4.0000,  4.4000,  5.0700])
RO_GAS = np.array([1e-4,    0.6000,  0.8000,  0.9165,  0.9798,  1.0000,
                   1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  0.9798,
                   0.9165,  0.8000,  0.6299,  0.4842,  0.3660,  0.2711,
                   0.1967,  0.1410,  0.1026,  0.0750,  0.0750])

Z_AL   = np.array([-3.5000, -3.4000, -3.3000, -3.1000, -2.9000, -2.7000,
                   -2.5000, -2.3000, -2.1000, -2.0000, -1.5000, -0.5000,
                    0.5000,  1.5000,  2.0000,  2.1000,  2.3000,  2.5000,
                    2.7000,  2.9000,  3.1000,  3.3000,  3.5000,  3.7000,
                    3.9000,  4.0000,  4.5000,  5.0000,  5.1000])
RO_AL  = np.array([0.0000,  0.3803,  0.5287,  0.7206,  0.8480,  0.9375,
                   0.9994,  1.0386,  1.0600,  1.0600,  1.0600,  1.0600,
                   1.0600,  1.0600,  1.0600,  1.0600,  1.0386,  0.9994,
                   0.9375,  0.8480,  0.7206,  0.5747,  0.4708,  0.4015,
                   0.3621,  0.3500,  0.3500,  0.3500,  0.3500])
RO_CFRP = np.where(RO_AL > 0, RO_AL + 0.09, 0.0)   # 0.9 mm wrap

he3_r       = 1.0    # gas bore radius [cm]
he3_total_r = 1.15   # CFRP outer radius at the barrel [cm]
he3_half_y  = (Z_AL[-1] - Z_AL[0]) / 2  # for plot extents


# ─────────────────────────────────────────────────────────────────────────────
# Arm geometry
# ─────────────────────────────────────────────────────────────────────────────
DIST_X   = CFG['mm_distance_x_cm']      # ±X arms (D, B) front-face distance
DIST_Z   = CFG['mm_distance_z_cm']      # ±Z arms (A, C) front-face distance
PINWHEEL = CFG['mm_pinwheel_shift_cm']  # per arm 0..3 = D,B,A,C; tangential along −u_hat
dist     = max(DIST_X, DIST_Z)          # representative extent (world sizing / limits)


def _arm(idx, label, w_hat, u_hat, d):
    """Arm def with two front-face centres:
       ff        — MM-shifted (pinwheel along −u_hat = +t)  [MM/plastics/LS]
       ff_struct — un-shifted mechanical structure centre    [SiPM wall]"""
    w = np.array(w_hat, dtype=float)
    u = np.array(u_hat, dtype=float)
    ff_struct = w * d
    ff = ff_struct + (-u) * PINWHEEL[idx]
    return dict(id=idx, label=label, ff=ff, ff_struct=ff_struct, u_hat=u, w_hat=w)


ARM_DEF = [
    _arm(0, '+X (D)', ( 1, 0, 0), (0, 0, -1), DIST_X),
    _arm(1, '−X (B)', (-1, 0, 0), (0, 0,  1), DIST_X),
    _arm(2, '+Z (A)', ( 0, 0, 1), (1, 0,  0), DIST_Z),
    _arm(3, '−Z (C)', ( 0, 0, -1), (-1, 0, 0), DIST_Z),
]
V_HAT = np.array([0, 1, 0])

HW_U = dict(mm=CFG['mm_size_u_cm']/2, ls=CFG['ls_size_u_cm']/2, sipm=bw/2)
HW_V = dict(mm=CFG['mm_size_v_cm']/2, ls=CFG['ls_size_v_cm']/2,
            sipm=CFG['sipm_size_v_cm']/2)

# u-offset to each plastic bar's centre (two bars side-by-side with a gap)
bsc_u_offset = bscTape_hu + bsc_gap / 2.0


def arm_box_world(arm, w_f, w_b, hu, hv, u_offset=0.0, ff=None):
    """World-frame centre and half-sizes for a slab in arm local coords.
    u_offset shifts the centre along the arm's u axis; ff overrides the front
    face (default = MM-shifted arm['ff']; pass arm['ff_struct'] for the SiPM)."""
    base = arm['ff'] if ff is None else ff
    wh, uh = arm['w_hat'], arm['u_hat']
    w_cen = (w_f + w_b) / 2
    hw    = (w_b - w_f) / 2
    centre = base + wh * w_cen + uh * u_offset
    hx = abs(wh[0])*hw + abs(uh[0])*hu + abs(V_HAT[0])*hv
    hy = abs(wh[1])*hw + abs(uh[1])*hu + abs(V_HAT[1])*hv
    hz = abs(wh[2])*hw + abs(uh[2])*hu + abs(V_HAT[2])*hv
    return centre, np.array([hx, hy, hz])


def _sipm_bar_offsets():
    """(index, u-offset) for all 20 bars; centered on the structure."""
    return [(i, bw * (i - (CFG['sipm_n_bars'] - 1) / 2.0))
            for i in range(CFG['sipm_n_bars'])]


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — Top-down 2D (matplotlib)
# ─────────────────────────────────────────────────────────────────────────────
def _rect_xz(ax, cen, hs, col, alpha, ls='-', lw=0.5, ec='k'):
    # screen mapping (adopted convention): x = Z (East →), y = X (North ↑)
    sx0, sy0 = cen[2] - hs[2], cen[0] - hs[0]
    ax.add_patch(Rectangle((sx0, sy0), 2*hs[2], 2*hs[0], linewidth=lw,
                           edgecolor=ec, facecolor=col, alpha=alpha,
                           linestyle=ls, zorder=2))


def plot_2d_topdown():
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')

    # He-3 capsule cross-section at the barrel — largest first, gas on top
    for r, fc, zo in [(1.15, '#404040', 3), (1.06, '#b0b0b0', 4), (1.00, '#99d8f5', 5)]:
        ax.add_patch(plt.Circle((0, 0), r, color=fc, zorder=zo))

    for arm in ARM_DEF:
        # MM + PCB (centered on the MM)
        for lname, w_f, w_b, col in MM_PCB_2D:
            cen, hs = arm_box_world(arm, w_f, w_b, HW_U['mm'], 0)
            _rect_xz(ax, cen, hs, col, 0.75)

        # SiPM wall — light container outline (3.3 cm, mostly empty) + thin bars
        # (3 mm scint) centered inside, on the STRUCTURE; 4 un-read → transparent
        cen, hs = arm_box_world(arm, w_sipm_f, w_sipm_b, sipm_wall_hu, 0,
                                ff=arm['ff_struct'])
        _rect_xz(ax, cen, hs, 'none', 1.0, lw=0.8, ec='0.6')
        for i, u_i in _sipm_bar_offsets():
            cen, hs = arm_box_world(arm, w_sipm_sc_f, w_sipm_sc_b, bw/2, 0,
                                    u_offset=u_i, ff=arm['ff_struct'])
            if i in SIPM_READOUT:
                _rect_xz(ax, cen, hs, C_SC, 0.90, lw=0.4)
            else:
                _rect_xz(ax, cen, hs, C_SC, 0.15, ls='--', lw=0.6, ec='0.6')

        # Plastics — two bars, centered on the MM
        for u_sign in (-1.0, +1.0):
            cen, hs = arm_box_world(arm, w_bsc_f, w_bsc_b, bscTape_hu, 0,
                                    u_offset=u_sign * bsc_u_offset)
            _rect_xz(ax, cen, hs, C_BS, 0.80)

        # LS box — single layer, centered on the MM
        for lname, w_f, w_b, col in LS_SUBS_2D:
            cen, hs = arm_box_world(arm, w_f, w_b, HW_U['ls'], 0)
            _rect_xz(ax, cen, hs, col, 0.75)

        # Arm label
        ff  = arm['ff']
        mid = ff + arm['w_hat'] * (stack_depth / 2)
        loff = ff / np.linalg.norm(ff) * 2.5
        ax.text(mid[2]+loff[2], mid[0]+loff[0],
                f"Arm {arm['id']}\n({arm['label']})",
                ha='center', va='center', fontsize=9, fontweight='bold', zorder=5,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='none', alpha=0.8))

    # Beam indicator (⊙ = +Y out of the page toward viewer)
    ax.annotate('beam ⊙\n(+Y, out of page)', xy=(0, 0), xytext=(2.5, 2.5),
                fontsize=8, zorder=6)

    cap_handles = [
        mpatches.Patch(color='#99d8f5', label='He-3 gas bore  (Ø20 × 40 mm + domed ends)'),
        mpatches.Patch(color='#404040', label='Capsule walls  (Al 0.6 mm + CFRP 0.9 mm)'),
    ]
    mm_handles = [
        mpatches.Patch(color=C_MM, alpha=0.75, label=f'MM drift gas  ({tDrift*10:.0f} mm drift)'),
        mpatches.Patch(color=C_PCB, alpha=0.75, label=f'PCB stack  ({t_PCB*10:.1f} mm)'),
    ]
    sc_handles = [
        mpatches.Patch(color=C_SC, alpha=0.85,
                       label=f'SiPM wall  ({CFG["sipm_n_readout"]}/{CFG["sipm_n_bars"]} bars '
                             f'× {bw:.1f} cm read out, on structure)'),
        mpatches.Patch(color=C_BS, alpha=0.80,
                       label=f'Plastics  (2×{int(bsc_u)}×{int(bsc_v)} cm, {bsc_th:.1f} cm, on MM)'),
    ]
    ls_handles = [
        mpatches.Patch(color=C_LS, alpha=0.75,
                       label=f'Liq. scint. LAB  (1×{tLS*10:.0f} mm, {CFG["ls_size_u_cm"]:.0f}×{CFG["ls_size_v_cm"]:.0f} cm, on MM)'),
        mpatches.Patch(color=C_CFRP, alpha=0.85, label='LS CFRP box walls + liners'),
    ]
    for handles, loc in [(cap_handles, 'upper left'), (mm_handles, 'upper right'),
                         (sc_handles, 'lower left'), (ls_handles, 'lower right')]:
        leg = ax.legend(handles=handles, loc=loc, fontsize=7.5, framealpha=0.85)
        ax.add_artist(leg)

    lim = dist + stack_depth + 4
    ax.set_xlim(-lim, lim);  ax.set_ylim(-lim, lim)
    ax.set_xlabel('Z  [cm]  (East →)', fontsize=12)
    ax.set_ylabel('X  [cm]  (North ↑)', fontsize=12)
    ax.set_title('MX17 Detector Geometry — Top-Down View\n'
                 'MM → SiPM wall → plastics → 1 LS layer   ·   '
                 '+Z right · +X up · beam +Y ⊙ out of page', fontsize=11.5)
    ax.axhline(0, color='0.7', lw=0.5, zorder=1)
    ax.axvline(0, color='0.7', lw=0.5, zorder=1)
    ax.grid(True, lw=0.3, alpha=0.5)

    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1b — Side view (XY plane): beam from below, capsule end-on
# ─────────────────────────────────────────────────────────────────────────────
def _profile_polygon(z, ro):
    return np.vstack([np.column_stack([ro, z]),
                      np.column_stack([-ro[::-1], z[::-1]])])


# Side-view layers for the ±X arms: (w_front, w_back, v_half, colour)
SIDE_LAYERS = [
    (w_MM_f,   w_MM_b,   HW_V['mm'],   C_MM),
    (w_PCB_f,  w_PCB_b,  HW_V['mm'],   C_PCB),
    (w_sipm_sc_f, w_sipm_sc_b, HW_V['sipm'], C_SC),   # thin 3 mm scint
    (w_bsc_f,  w_bsc_b,  bscTape_hv,   C_BS),
    (_lsDivFront_f, _lsDivFront_b, HW_V['ls'], C_CFRP),
    (_wLS1_f,       _wLS1_b,       HW_V['ls'], C_LS),
    (_lsDivRear_f,  _lsDivRear_b,  HW_V['ls'], C_CFRP),
]


def plot_2d_sideview():
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.set_aspect('equal')

    for z, ro, fc, zo, lab in [
        (Z_AL,  RO_CFRP, '#404040', 3, 'CFRP wrap (0.9 mm)'),
        (Z_AL,  RO_AL,   '#b0b0b0', 4, 'Al vessel (0.6 mm barrel; domed nose; Ø7 mm neck/valve)'),
        (Z_GAS, RO_GAS,  '#99d8f5', 5, 'He-3 gas (Ø20 mm × 40 mm + domed ends)'),
    ]:
        ax.add_patch(plt.Polygon(_profile_polygon(z, ro), closed=True,
                                 fc=fc, ec='k', lw=0.4, zorder=zo, label=lab))

    for sign in (+1, -1):
        for w_f, w_b, hv, col in SIDE_LAYERS:
            x0 = sign * (dist + w_f) if sign > 0 else sign * (dist + w_b)
            ax.add_patch(Rectangle((x0, -hv), w_b - w_f, 2 * hv,
                                   linewidth=0.4, edgecolor='k',
                                   facecolor=col, alpha=0.75, zorder=2))
        # SiPM container outline (3.3 cm, mostly empty)
        xf = sign * (dist + w_sipm_f) if sign > 0 else sign * (dist + w_sipm_b)
        ax.add_patch(Rectangle((xf, -HW_V['sipm']), w_sipm_b - w_sipm_f,
                               2 * HW_V['sipm'], fill=False, edgecolor='0.6',
                               lw=0.8, zorder=2))
        mid = sign * (dist + stack_depth / 2)
        ax.text(mid, HW_V['sipm'] + 3, f"Arm {0 if sign > 0 else 1} "
                f"({'+' if sign > 0 else '−'}X)", ha='center', fontsize=9,
                fontweight='bold')

    ax.annotate('', xy=(0, -5), xytext=(0, -20),
                arrowprops=dict(arrowstyle='-|>', color='firebrick', lw=2.5))
    ax.text(1.2, -14, 'beam (+Y)', color='firebrick', fontsize=10)

    lim = dist + stack_depth + 4
    ax.set_xlim(-lim, lim); ax.set_ylim(-32, 34)
    ax.set_xlabel('X  [cm]'); ax.set_ylabel('Y (beam)  [cm]')
    ax.set_title('MX17 Geometry — Side View (XY plane)\n'
                 'Beam from below hits the capsule end-on; arms 2,3 (±Z) not shown')
    ax.legend(loc='lower left', fontsize=8, framealpha=0.85)
    ax.grid(True, lw=0.3, alpha=0.5)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — 3D view (pyvista/VTK)
# ─────────────────────────────────────────────────────────────────────────────
def _add_box(pl, cen, hs, col, alpha, label=None):
    cx, cy, cz = cen
    hx, hy, hz = hs
    box = pv.Box(bounds=(cx-hx, cx+hx, cy-hy, cy+hy, cz-hz, cz+hz))
    pl.add_mesh(box, color=col, opacity=alpha, smooth_shading=True,
                show_edges=False, label=label)


def plot_3d_pyvista(out_path=None, interactive=True):
    pl = pv.Plotter(off_screen=not interactive, window_size=(1600, 1200))
    pl.set_background('white')
    pl.enable_depth_peeling(number_of_peels=8, occlusion_ratio=0.0)

    seen = set()

    def once(key, label):
        if key in seen:
            return None
        seen.add(key)
        return label

    for arm in ARM_DEF:
        # MM + PCB (on MM)
        cen, hs = arm_box_world(arm, w_MM_f, w_MM_b, HW_U['mm'], HW_V['mm'])
        _add_box(pl, cen, hs, (0.29, 0.56, 0.85), 0.80,
                 once('MM', f'MM drift gas  ({tDrift*10:.0f} mm drift)'))
        cen, hs = arm_box_world(arm, w_PCB_f, w_PCB_b, HW_U['mm'], HW_V['mm'])
        _add_box(pl, cen, hs, (0.36, 0.72, 0.36), 0.80,
                 once('PCB', f'PCB stack  ({t_PCB*10:.1f} mm)'))

        # SiPM wall — 16 read bars solid, 4 un-read faint (on structure)
        for i, u_i in _sipm_bar_offsets():
            cen, hs = arm_box_world(arm, w_sipm_sc_f, w_sipm_sc_b, bw/2,
                                    HW_V['sipm'], u_offset=u_i, ff=arm['ff_struct'])
            if i in SIPM_READOUT:
                _add_box(pl, cen, hs, (0.94, 0.75, 0.25), 0.85,
                         once('SC', f'SiPM wall  ({CFG["sipm_n_readout"]}/{CFG["sipm_n_bars"]} bars read out)'))
            else:
                _add_box(pl, cen, hs, (0.94, 0.75, 0.25), 0.12,
                         once('SCx', 'SiPM bar (not read out)'))
        # SiPM container outline (3.3 cm, mostly empty) — wireframe
        cen, hs = arm_box_world(arm, w_sipm_f, w_sipm_b, sipm_wall_hu,
                                HW_V['sipm'], ff=arm['ff_struct'])
        cx, cy, cz = cen
        hx, hy, hz = hs
        cbox = pv.Box(bounds=(cx-hx, cx+hx, cy-hy, cy+hy, cz-hz, cz+hz))
        pl.add_mesh(cbox, color=(0.5, 0.5, 0.5), opacity=0.3, style='wireframe',
                    line_width=1.5,
                    label=once('SCcont', 'SiPM container (3.3 cm, mostly empty)'))

        # Plastics — two bars (on MM)
        for u_sign in (-1.0, +1.0):
            cen, hs = arm_box_world(arm, w_bsc_f, w_bsc_b, bscTape_hu,
                                    bscTape_hv, u_offset=u_sign * bsc_u_offset)
            _add_box(pl, cen, hs, (0.88, 0.47, 0.13), 0.85,
                     once('BSC', f'Plastics  (2×{int(bsc_u)}×{int(bsc_v)} cm, {bsc_th:.1f} cm)'))

        # LS box (on MM)
        for lname, w_f, w_b, col in [
            ('LS_CFRP', _lsDivFront_f, _lsDivFront_b, (0.20, 0.20, 0.20)),
            ('LS',      _wLS1_f,       _wLS1_b,       (0.85, 0.33, 0.31)),
            ('LS_CFRP', _lsDivRear_f,  _lsDivRear_b,  (0.20, 0.20, 0.20)),
        ]:
            cen, hs = arm_box_world(arm, w_f, w_b, HW_U['ls'], HW_V['ls'])
            lbl = None
            if lname == 'LS':
                lbl = once('LS', f'Liq. scint. LAB  (1×{tLS*10:.0f} mm, {CFG["ls_size_u_cm"]:.0f}×{CFG["ls_size_v_cm"]:.0f} cm)')
            elif lname == 'LS_CFRP':
                lbl = once('LSC', 'LS CFRP box walls + liners')
            _add_box(pl, cen, hs, col, 0.85 if lname == 'LS_CFRP' else 0.70, lbl)

    # He-3 capsule
    def polycone_mesh(z, ro):
        pts = np.column_stack([ro, np.zeros_like(ro), z])
        line = pv.lines_from_points(pts)
        mesh = line.extrude_rotate(resolution=72, capping=True)
        mesh.rotate_x(-90.0, inplace=True)
        return mesh

    caps = [
        (Z_AL,  RO_CFRP, (0.16, 0.16, 0.16), 0.55, 'CFRP wrap (0.9 mm)'),
        (Z_AL,  RO_AL,   (0.67, 0.67, 0.67), 0.65,
         'Al vessel (0.6 mm barrel; domed nose; Ø7 mm neck/valve)'),
        (Z_GAS, RO_GAS,  (0.60, 0.85, 0.96), 0.90,
         'He-3 gas (Ø20 mm × 40 mm + domed ends)'),
    ]
    for z, ro, col, alpha, llabel in caps:
        pl.add_mesh(polycone_mesh(z, ro), color=col, opacity=alpha,
                    smooth_shading=True, show_edges=False, label=llabel)

    yl = he3_half_y + 10
    arrow = pv.Arrow(start=(0, -yl, 0), direction=(0, 1, 0),
                     scale=2*yl, tip_length=0.08, tip_radius=0.025, shaft_radius=0.008)
    pl.add_mesh(arrow, color='firebrick', label='+Y beam axis')
    pl.add_point_labels([(0, yl + 1, 0)], ['+Y (beam)'], text_color='firebrick',
                        font_size=14, bold=True, show_points=False, always_visible=True)

    for arm in ARM_DEF:
        pt = arm['ff'] + arm['w_hat'] * (stack_depth + 4)
        pl.add_point_labels([tuple(pt)], [f"Arm {arm['id']} ({arm['label']})"],
                            text_color='black', font_size=11, bold=True,
                            show_points=False, always_visible=True)

    pl.add_legend(size=(0.30, 0.24), loc='upper left',
                  bcolor=(0.95, 0.95, 0.95), border=True, background_opacity=0.9)
    pl.add_axes(xlabel='X [cm]', ylabel='Y [cm]', zlabel='Z [cm]',
                line_width=3, labels_off=False)
    pl.add_title('MX17 Detector Geometry — 3D View\n'
                 'MM → SiPM wall → plastics → 1 LS layer  |  Beam along +Y',
                 font_size=12, color='black')

    lim = dist + stack_depth + 5
    cam_az_deg = 227.0
    cam_el     = 1.55
    r_xz       = 3.96 * lim
    az = np.radians(cam_az_deg)
    pl.camera_position = [
        (r_xz*np.cos(az), lim*cam_el, r_xz*np.sin(az)),
        (0, 0, 0),
        (0, 1, 0),
    ]

    if interactive:
        pl.show()
    else:
        pl.screenshot(out_path, return_img=False)
        pl.close()
        print(f"Saved: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    import sys
    interactive = '--interactive' in sys.argv
    no_3d       = '--no-3d' in sys.argv or pv is None
    _here = os.path.dirname(os.path.abspath(__file__))

    print(f"MM front face distance : ±X (D,B) {DIST_X:.2f} cm | ±Z (A,C) {DIST_Z:.2f} cm")
    print(f"Pinwheel shift [D,B,A,C]: {PINWHEEL} cm")
    print(f"MM+PCB back            : {(w_PCB_b):.2f} cm")
    print(f"SiPM container         : {w_sipm_f:.2f} → {w_sipm_b:.2f} cm  "
          f"(read out bars {sorted(SIPM_READOUT)[0]}–{sorted(SIPM_READOUT)[-1]} of "
          f"0–{CFG['sipm_n_bars']-1})")
    print(f"Plastics               : {w_bsc_f:.2f} → {w_bsc_b:.2f} cm  ({bsc_th:.1f} cm thick)")
    print(f"LS box                 : {w_LS_f:.2f} → {w_LS_b:.2f} cm  (single {tLS*10:.0f} mm LAB layer)")
    print(f"Total stack depth      : {stack_depth:.2f} cm")
    print(f"Arm outer edge         : {dist + stack_depth:.1f} cm from origin")

    out1 = os.path.join(_here, 'mx17_geometry_topdown.png')
    out1b = os.path.join(_here, 'mx17_geometry_sideview')
    out2 = os.path.join(_here, 'mx17_geometry_3d.png')

    fig1 = plot_2d_topdown()
    fig1.savefig(out1, dpi=150, bbox_inches='tight')
    print(f"Saved: {out1}")

    fig1b = plot_2d_sideview()
    fig1b.savefig(out1b + '.png', dpi=150, bbox_inches='tight')
    fig1b.savefig(out1b + '.pdf', bbox_inches='tight')
    print(f"Saved: {out1b}.png/.pdf")

    if no_3d:
        print("(3D view skipped: --no-3d or pyvista unavailable)")
    else:
        plot_3d_pyvista(out_path=out2, interactive=interactive)
    if interactive:
        plt.show()
