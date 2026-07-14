#!/usr/bin/env python3
"""
plot_geometry.py
Visualise the MX17 4-arm detector geometry.
  Figure 1 — top-down 2D (+Z right, +X up, beam +Y ⊙ out of page)  [matplotlib]
  Figure 2 — 3D isometric view                                   [pyvista/VTK]

No Geant4 dependency; dimensions are taken directly from SimConfig defaults.
Run from anywhere:  python scripts/plot_geometry.py
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
    scint_size_u_cm       = 48.0,   # trigger scint wall
    scint_size_v_cm       = 48.0,
    ls_size_u_cm          = 45.0,   # liquid scintillator (2 layers, 45×45 cm)
    ls_size_v_cm          = 45.0,
    ls_cfrp_mm            =  2.0,   # structural CFRP wall thickness
    ls_inner_cfrp_um      = 600.0,  # inner CFRP liner before each LAB layer [µm]
    ls_inner_al_um        =  40.0,  # Al liner before each LAB layer [µm]
    ls_thick_cm           =  2.0,   # LAB layer thickness
    backscint_u_cm        = 20.0,   # each bar (two per arm, side-by-side)
    backscint_v_cm        = 30.0,
    backscint_thick_cm    =  2.0,
    backscint_gap_cm      =  0.3,   # gap between wrapped bars
    backscint_tape_um     = 200.0,  # black mylar tape (outermost)
    backscint_al_um       =  20.0,  # Al foil on scintillator surface [µm]
    gap_pcb_to_scint_mm   = 20.0,
    gap_scint_to_ls_mm    = 20.0,
    gap_ls_to_back_mm     =  1.0,   # gap_ls_to_backscint_mm in SimConfig.hh
)

# ─────────────────────────────────────────────────────────────────────────────
# Layer thicknesses  [cm]
# ─────────────────────────────────────────────────────────────────────────────
_um = 1e-4
_mm = 0.1

tMylar    =  40 * _um
tAlWin    = 0.1 * _um
tKapCath  =  50 * _um
tCuCath   =   9 * _um
tDrift    =  30 * _mm
tMesh     =  30 * _um
tAmp      = 150 * _um
tResPaste = 100 * _um

tPCBKap  =  50 * _um
tPCBCu   =  26 * _um
tPCBFR4  = 100 * _um
tPCBRoh  =   5 * _mm
tPCBAl   =  50 * _um

tBlkTape = 200 * _um    # black mylar tape (trigger scint wall + back scint wrapping)
tPlScint =   3 * _mm   # trigger scint wall
tScAl    =  50 * _um

tLSCfrp      = CFG['ls_cfrp_mm']        * _mm   # 2 mm structural CFRP wall
tLSInnerCfrp = CFG['ls_inner_cfrp_um'] * _um   # 600 µm inner CFRP liner
tLSInnerAl   = CFG['ls_inner_al_um']   * _um   # 40 µm Al liner
tLS          = CFG['ls_thick_cm']      * 1.0   # 2 cm LAB layer (already in cm)
tTape        = CFG['backscint_tape_um'] * _um   # back scint tape (outer)
tBscAl       = CFG['backscint_al_um']  * _um   # back scint Al foil

t_MM     = tMylar + tAlWin + tKapCath + tCuCath + tDrift + tMesh + tAmp + tResPaste
t_PCB    = tPCBKap + 4*(tPCBCu + tPCBFR4) + tPCBRoh + tPCBAl
t_scint  = 2*tBlkTape + tPlScint + tScAl
t_LS     = 3*tLSCfrp + 2*(tLSInnerCfrp + tLSInnerAl + tLS)   # full LS stack

gap2     = CFG['gap_pcb_to_scint_mm'] * _mm
gap3     = CFG['gap_scint_to_ls_mm']  * _mm
gap4     = CFG['gap_ls_to_back_mm']   * _mm

# Back scint bar: tape (outer) → Al foil → PVT
bsc_u   = CFG['backscint_u_cm']
bsc_v   = CFG['backscint_v_cm']
bsc_th  = CFG['backscint_thick_cm']
bsc_gap = CFG['backscint_gap_cm']
bscTape_hu = (bsc_u  + 2*tBscAl + 2*tTape) / 2   # outermost envelope half-width
bscTape_hv = (bsc_v  + 2*tBscAl + 2*tTape) / 2
bscTape_hw = (bsc_th + 2*tBscAl + 2*tTape) / 2
t_bsc    = 2 * bscTape_hw

w_MM_f   = 0.0;               w_MM_b   = t_MM
w_PCB_f  = w_MM_b;            w_PCB_b  = w_PCB_f  + t_PCB
w_sc_f   = w_PCB_b + gap2;    w_sc_b   = w_sc_f   + t_scint
w_LS_f   = w_sc_b  + gap3;    w_LS_b   = w_LS_f   + t_LS
w_bsc_f  = w_LS_b  + gap4;    w_bsc_b  = w_bsc_f  + t_bsc
stack_depth = w_bsc_b

# LS sub-layer positions (within the LS stack):
#   CFRP_1 | InnerCFRP_1+Al_1 | LAB_1 | CFRP_2 | InnerCFRP_2+Al_2 | LAB_2 | CFRP_3
_lsDivFront_f = w_LS_f
_lsDivFront_b = w_LS_f + tLSCfrp + tLSInnerCfrp + tLSInnerAl   # front wall + liners
_wLS1_f       = _lsDivFront_b
_wLS1_b       = _wLS1_f + tLS
_lsDivMid_f   = _wLS1_b
_lsDivMid_b   = _wLS1_b + tLSCfrp + tLSInnerCfrp + tLSInnerAl  # middle wall + liners
_wLS2_f       = _lsDivMid_b
_wLS2_b       = _wLS2_f + tLS
_lsDivRear_f  = _wLS2_b
_lsDivRear_b  = w_LS_b                                            # rear CFRP wall


# ─────────────────────────────────────────────────────────────────────────────
# Arm geometry
# ─────────────────────────────────────────────────────────────────────────────
DIST_X   = CFG['mm_distance_x_cm']      # ±X arms (D, B) front-face distance
DIST_Z   = CFG['mm_distance_z_cm']      # ±Z arms (A, C) front-face distance
PINWHEEL = CFG['mm_pinwheel_shift_cm']  # per arm 0..3 = D,B,A,C; tangential along −u_hat
dist     = max(DIST_X, DIST_Z)          # representative extent (world sizing / limits)


def _arm(idx, label, w_hat, u_hat, d):
    """Arm def with front-face centre shifted tangentially (pinwheel, −u_hat)."""
    w = np.array(w_hat, dtype=float)
    u = np.array(u_hat, dtype=float)
    ff = w * d + (-u) * PINWHEEL[idx]
    return dict(id=idx, label=label, ff=ff, u_hat=u, w_hat=w)


ARM_DEF = [
    _arm(0, '+X (D)', ( 1, 0, 0), (0, 0, -1), DIST_X),
    _arm(1, '−X (B)', (-1, 0, 0), (0, 0,  1), DIST_X),
    _arm(2, '+Z (A)', ( 0, 0, 1), (1, 0,  0), DIST_Z),
    _arm(3, '−Z (C)', ( 0, 0, -1), (-1, 0, 0), DIST_Z),
]
V_HAT = np.array([0, 1, 0])

HW_U = dict(mm=CFG['mm_size_u_cm']/2, scint=CFG['scint_size_u_cm']/2,
            ls=CFG['ls_size_u_cm']/2)
HW_V = dict(mm=CFG['mm_size_v_cm']/2, scint=CFG['scint_size_v_cm']/2,
            ls=CFG['ls_size_v_cm']/2)

# Back scint bars are drawn individually (not as a combined bounding box)
# so we can show the gap. This constant is the u-offset to each bar's centre.
bsc_u_offset = bscTape_hu + bsc_gap / 2.0   # 12.52 + 0.15 = 12.67 cm

# Standard slab layers for 2D/3D
# Back scints are omitted here and handled separately in each plot function.
LAYERS_2D = [
    ('MM',        w_MM_f,        w_MM_b,        'mm',    'mm',    '#4a90d9'),
    ('PCB',       w_PCB_f,       w_PCB_b,       'mm',    'mm',    '#5cb85c'),
    ('PlScint',   w_sc_f,        w_sc_b,        'scint', 'scint', '#f0c040'),
    ('LS_CFRP',   _lsDivFront_f, _lsDivFront_b, 'ls',    'ls',    '#303030'),
    ('LiqScint1', _wLS1_f,       _wLS1_b,       'ls',    'ls',    '#d9534f'),
    ('LS_CFRP',   _lsDivMid_f,   _lsDivMid_b,   'ls',    'ls',    '#303030'),
    ('LiqScint2', _wLS2_f,       _wLS2_b,       'ls',    'ls',    '#d9534f'),
    ('LS_CFRP',   _lsDivRear_f,  _lsDivRear_b,  'ls',    'ls',    '#303030'),
]

# ─────────────────────────────────────────────────────────────────────────────
# He-3 capsule — full STEP-derived polycone profiles [cm]
# (axis along beam = Y; keep in sync with src/DetectorConstruction.cc)
#
# Profiles extracted from the CAD master
#   "MASTINU X17 HPRV 00 01 (Cylinder D20 L40 mm).step"
# by sectioning the (axisymmetric) Al solid at successive Y levels and
# recording the outer wall radius and the inner gas-cavity radius.
# Real dimensions recovered from the STEP:
#   • He-3 gas bore : r = 10.0 mm, 40 mm cylinder (Y = ±20 mm)
#                     + lower hemispherical end + upper conical fill channel
#                     necking to the Ø1.5 mm valve bore (r = 0.75 mm)
#   • Al vessel     : 0.6 mm barrel wall (OD 21.2 mm); domed solid nose
#                     (Y −35 → −20 mm); shoulder taper into the Ø7 mm
#                     neck/valve (r = 3.5 mm), full length Y −35 → +51 mm
#   • CFRP overwrap : 0.9 mm radial, applied over the whole vessel
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


def arm_box_world(arm, w_f, w_b, hu, hv, u_offset=0.0):
    """World-frame centre and half-sizes for a slab in arm local coords.
    u_offset shifts the centre along the arm's u axis (default 0 = arm centre)."""
    ff, wh, uh = arm['ff'], arm['w_hat'], arm['u_hat']
    w_cen = (w_f + w_b) / 2
    hw    = (w_b - w_f) / 2
    centre = ff + wh * w_cen + uh * u_offset
    hx = abs(wh[0])*hw + abs(uh[0])*hu + abs(V_HAT[0])*hv
    hy = abs(wh[1])*hw + abs(uh[1])*hu + abs(V_HAT[1])*hv
    hz = abs(wh[2])*hw + abs(uh[2])*hu + abs(V_HAT[2])*hv
    return centre, np.array([hx, hy, hz])


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — Top-down 2D (matplotlib)
# ─────────────────────────────────────────────────────────────────────────────
def plot_2d_topdown():
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')

    # He-3 capsule cross-section at the barrel — largest first, gas on top
    for r, fc, zo in [
        (1.15, '#404040', 3),   # CFRP outer
        (1.06, '#b0b0b0', 4),   # Al outer
        (1.00, '#99d8f5', 5),   # gas bore
    ]:
        ax.add_patch(plt.Circle((0, 0), r, color=fc, zorder=zo))

    # Arm detector slabs projected onto XZ
    for arm in ARM_DEF:
        for lname, w_f, w_b, ukey, vkey, col in LAYERS_2D:
            cen, hs = arm_box_world(arm, w_f, w_b, HW_U[ukey], 0)
            # screen mapping (adopted convention): x = Z (East →), y = X (North ↑)
            sx0, sy0 = cen[2] - hs[2], cen[0] - hs[0]
            ax.add_patch(Rectangle((sx0, sy0), 2*hs[2], 2*hs[0],
                                    linewidth=0.5, edgecolor='k',
                                    facecolor=col, alpha=0.75, zorder=2))
        ff  = arm['ff']
        mid = ff + arm['w_hat'] * (stack_depth / 2)
        loff = ff / np.linalg.norm(ff) * 2.5
        ax.text(mid[2]+loff[2], mid[0]+loff[0],
                f"Arm {arm['id']}\n({arm['label']})",
                ha='center', va='center', fontsize=9, fontweight='bold', zorder=5,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='none', alpha=0.8))

    # Back plastic scints — two individual bars per arm with gap visible
    for arm in ARM_DEF:
        for u_sign in [-1.0, +1.0]:
            cen, hs = arm_box_world(arm, w_bsc_f, w_bsc_b,
                                    bscTape_hu, bscTape_hv,
                                    u_offset=u_sign * bsc_u_offset)
            sx0, sy0 = cen[2] - hs[2], cen[0] - hs[0]   # screen: x=Z, y=X
            ax.add_patch(Rectangle((sx0, sy0), 2*hs[2], 2*hs[0],
                                    linewidth=0.5, edgecolor='k',
                                    facecolor='#e07820', alpha=0.80, zorder=2))

    # Beam indicator (⊙ = +Y out of the page toward viewer)
    ax.annotate('beam ⊙\n(+Y, out of page)', xy=(0, 0), xytext=(2.5, 2.5),
                fontsize=8, zorder=6)

    # Legend split into four corner boxes (grouped by subsystem) so it sits in
    # the empty diagonal corners rather than overlapping the cardinal arms.
    cap_handles = [
        mpatches.Patch(color='#99d8f5',
                       label='He-3 gas bore  (Ø20 × 40 mm + domed ends)'),
        mpatches.Patch(color='#404040',
                       label='Capsule walls  (Al 0.6 mm + CFRP 0.9 mm)'),
    ]
    mm_handles = [
        mpatches.Patch(color='#4a90d9', alpha=0.75,
                       label=f'MM drift gas  ({tDrift*10:.0f} mm drift)'),
        mpatches.Patch(color='#5cb85c', alpha=0.75,
                       label=f'PCB stack  ({t_PCB*10:.1f} mm)'),
    ]
    sc_handles = [
        mpatches.Patch(color='#f0c040', alpha=0.75,
                       label=f'Trigger plastic scint  ({tPlScint*10:.0f} mm, 48×48 cm)'),
        mpatches.Patch(color='#e07820', alpha=0.75,
                       label=f'Back plastic scints  ({int(bsc_u)} mm, {int(bsc_v)}×(2×{int(bsc_th*10)}) cm)'),
    ]
    ls_handles = [
        mpatches.Patch(color='#d9534f', alpha=0.75,
                       label=f'Liq. scint. LAB  (2×{tLS*10:.0f} mm, {CFG["ls_size_u_cm"]:.0f}×{CFG["ls_size_v_cm"]:.0f} cm)'),
        mpatches.Patch(color='#303030', alpha=0.85,
                       label=f'LS CFRP walls + liners  ({tLSCfrp*10:.0f} mm + {tLSInnerCfrp*1e3:.0f} µm + {tLSInnerAl*1e3:.0f} µm Al)'),
    ]
    for handles, loc in [(cap_handles, 'upper left'),
                         (mm_handles,  'upper right'),
                         (sc_handles,  'lower left'),
                         (ls_handles,  'lower right')]:
        leg = ax.legend(handles=handles, loc=loc, fontsize=7.5, framealpha=0.85)
        ax.add_artist(leg)

    lim = dist + stack_depth + 4
    ax.set_xlim(-lim, lim);  ax.set_ylim(-lim, lim)
    ax.set_xlabel('Z  [cm]  (East →)', fontsize=12)
    ax.set_ylabel('X  [cm]  (North ↑)', fontsize=12)
    ax.set_title('MX17 Detector Geometry — Top-Down View\n'
                 '+Z right · +X up · beam +Y ⊙ out of page (right-handed)', fontsize=12)
    ax.axhline(0, color='0.7', lw=0.5, zorder=1)
    ax.axvline(0, color='0.7', lw=0.5, zorder=1)
    ax.grid(True, lw=0.3, alpha=0.5)

    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1b — Side view (XY plane): beam from below, capsule end-on
# ─────────────────────────────────────────────────────────────────────────────
def _profile_polygon(z, ro):
    """Closed polygon (x=±r, y=z) for a polycone profile (axis along Y)."""
    return np.vstack([np.column_stack([ro, z]),
                      np.column_stack([-ro[::-1], z[::-1]])])


def plot_2d_sideview():
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.set_aspect('equal')

    # Capsule: CFRP → Al → gas (axis along Y = beam)
    for z, ro, fc, zo, lab in [
        (Z_AL,  RO_CFRP, '#404040', 3, 'CFRP wrap (0.9 mm)'),
        (Z_AL,  RO_AL,   '#b0b0b0', 4, 'Al vessel (0.6 mm barrel; '
                                       'domed nose; Ø7 mm neck/valve)'),
        (Z_GAS, RO_GAS,  '#99d8f5', 5, 'He-3 gas (Ø20 mm × 40 mm + domed ends)'),
    ]:
        ax.add_patch(plt.Polygon(_profile_polygon(z, ro), closed=True,
                                 fc=fc, ec='k', lw=0.4, zorder=zo, label=lab))

    # Arms 0 (+X) and 1 (−X) in side projection (v extent along Y)
    for sign in (+1, -1):
        for lname, w_f, w_b, ukey, vkey, col in LAYERS_2D:
            hv = HW_V[vkey]
            x0 = sign * (dist + w_f) if sign > 0 else sign * (dist + w_b)
            ax.add_patch(Rectangle((x0, -hv), w_b - w_f, 2 * hv,
                                   linewidth=0.4, edgecolor='k',
                                   facecolor=col, alpha=0.75, zorder=2))
        for u_sign in (-1.0, +1.0):   # back scint bars (v = 30 cm)
            x0 = sign * (dist + w_bsc_f) if sign > 0 else sign * (dist + w_bsc_b)
            ax.add_patch(Rectangle((x0, -bscTape_hv), t_bsc, 2 * bscTape_hv,
                                   linewidth=0.4, edgecolor='k',
                                   facecolor='#e07820', alpha=0.8, zorder=2))
        mid = sign * (dist + stack_depth / 2)
        ax.text(mid, HW_V['scint'] + 3, f"Arm {0 if sign > 0 else 1} "
                f"({'+' if sign > 0 else '−'}X)", ha='center', fontsize=9,
                fontweight='bold')

    # Beam arrow from below along +Y
    ax.annotate('', xy=(0, -5), xytext=(0, -20),
                arrowprops=dict(arrowstyle='-|>', color='firebrick', lw=2.5))
    ax.text(1.2, -14, 'beam (+Y)', color='firebrick', fontsize=10)

    lim = dist + stack_depth + 4
    ax.set_xlim(-lim, lim); ax.set_ylim(-28, 30)
    ax.set_xlabel('X  [cm]'); ax.set_ylabel('Y (beam)  [cm]')
    ax.set_title('MX17 Geometry — Side View (XY plane)\n'
                 'Beam from below hits the capsule end-on; '
                 'arms 2,3 (±Z) not shown')
    ax.legend(loc='lower left', fontsize=8, framealpha=0.85)
    ax.grid(True, lw=0.3, alpha=0.5)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — 3D view (pyvista/VTK)
# ─────────────────────────────────────────────────────────────────────────────
def plot_3d_pyvista(out_path=None, interactive=True):
    pl = pv.Plotter(off_screen=not interactive, window_size=(1600, 1200))
    pl.set_background('white')

    # Depth peeling gives correct transparency sorting for overlapping meshes
    pl.enable_depth_peeling(number_of_peels=8, occlusion_ratio=0.0)

    # Layer definitions: (name, w_front, w_back, hw_u_key, RGB, opacity, legend_label)
    layer_info = [
        ('MM',  w_MM_f,  w_MM_b,  'mm',    (0.29, 0.56, 0.85), 0.80,
         f'MM drift gas  ({tDrift*10:.0f} mm drift)'),
        ('PCB', w_PCB_f, w_PCB_b, 'mm',    (0.36, 0.72, 0.36), 0.80,
         f'PCB stack  ({t_PCB*10:.1f} mm)'),
        ('SC',  w_sc_f,  w_sc_b,  'scint', (0.94, 0.75, 0.25), 0.85,
         f'Trigger plastic scint  ({tPlScint*10:.0f} mm, 48×48 cm)'),
        ('LS_CFRP', _lsDivFront_f, _lsDivFront_b, 'ls', (0.20, 0.20, 0.20), 0.85,
         f'LS CFRP walls + liners  ({tLSCfrp*10:.0f} mm + {tLSInnerCfrp*1e3:.0f} µm + {tLSInnerAl*1e3:.0f} µm Al)'),
        ('LS',  _wLS1_f, _wLS1_b, 'ls', (0.85, 0.33, 0.31), 0.70,
         f'Liq. scint. LAB layers  (2×{tLS*10:.0f} mm, {CFG["ls_size_u_cm"]:.0f}×{CFG["ls_size_v_cm"]:.0f} cm)'),
        ('LS_CFRP', _lsDivMid_f,  _lsDivMid_b,  'ls', (0.20, 0.20, 0.20), 0.85, None),
        ('LS',  _wLS2_f, _wLS2_b, 'ls', (0.85, 0.33, 0.31), 0.70, None),
        ('LS_CFRP', _lsDivRear_f, _lsDivRear_b, 'ls', (0.20, 0.20, 0.20), 0.85, None),
    ]
    bsc_col   = (0.88, 0.47, 0.13)
    bsc_label = f'Back plastic scints  ({int(bsc_u)}mm, {int(bsc_v)}×(2×{int(bsc_th*10)}) cm)'

    seen = set()
    for arm in ARM_DEF:
        # Standard full-width slabs
        for lname, w_f, w_b, ukey, col, alpha, llabel in layer_info:
            cen, hs = arm_box_world(arm, w_f, w_b, HW_U[ukey], HW_V[ukey])
            cx, cy, cz = cen
            hx, hy, hz = hs
            box = pv.Box(bounds=(cx-hx, cx+hx, cy-hy, cy+hy, cz-hz, cz+hz))
            pl.add_mesh(box, color=col, opacity=alpha,
                        smooth_shading=True, show_edges=False,
                        label=llabel if lname not in seen else None)
            seen.add(lname)

        # Back scint bars — two per arm at ±u_offset
        for u_sign in [-1.0, +1.0]:
            cen, hs = arm_box_world(arm, w_bsc_f, w_bsc_b,
                                    bscTape_hu, bscTape_hv,
                                    u_offset=u_sign * bsc_u_offset)
            cx, cy, cz = cen
            hx, hy, hz = hs
            box = pv.Box(bounds=(cx-hx, cx+hx, cy-hy, cy+hy, cz-hz, cz+hz))
            pl.add_mesh(box, color=bsc_col, opacity=0.85,
                        smooth_shading=True, show_edges=False,
                        label=bsc_label if 'BSC' not in seen else None)
            seen.add('BSC')

    # He-3 capsule: STEP polycone profiles revolved about the beam (Y) axis
    def polycone_mesh(z, ro):
        pts = np.column_stack([ro, np.zeros_like(ro), z])
        line = pv.lines_from_points(pts)
        mesh = line.extrude_rotate(resolution=72, capping=True)
        mesh.rotate_x(-90.0, inplace=True)   # local Z (axis) → world Y (beam)
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

    # Beam axis arrow along +Y
    yl = he3_half_y + 10
    arrow = pv.Arrow(start=(0, -yl, 0), direction=(0, 1, 0),
                     scale=2*yl, tip_length=0.08, tip_radius=0.025,
                     shaft_radius=0.008)
    pl.add_mesh(arrow, color='firebrick', label='+Y beam axis')
    pl.add_point_labels(
        [(0, yl + 1, 0)], ['+Y (beam)'],
        text_color='firebrick', font_size=14, bold=True,
        show_points=False, always_visible=True,
    )

    # Arm labels
    for arm in ARM_DEF:
        pt = arm['ff'] + arm['w_hat'] * (stack_depth + 4)
        pl.add_point_labels(
            [tuple(pt)], [f"Arm {arm['id']} ({arm['label']})"],
            text_color='black', font_size=11, bold=True,
            show_points=False, always_visible=True,
        )

    # Legend
    pl.add_legend(size=(0.28, 0.22), loc='upper left',
                  bcolor=(0.95, 0.95, 0.95), border=True,
                  background_opacity=0.9)

    # Axes indicator (bottom-left corner)
    pl.add_axes(xlabel='X [cm]', ylabel='Y [cm]', zlabel='Z [cm]',
                line_width=3, labels_off=False)

    # Title
    pl.add_title('MX17 Detector Geometry — 3D View\n'
                 'Arms 0,1 at ±X  |  Arms 2,3 at ±Z  |  Beam along +Y',
                 font_size=12, color='black')

    # Camera: look down a diagonal gap between two arms at the target; Y is "up".
    # Azimuth 225° = straight down the −X,−Z gap; rotated a little to the left so
    # the He-3 target sits centred in the gap between the detectors.
    lim = dist + stack_depth + 5
    cam_az_deg = 227.0            # XZ-plane azimuth of the camera (225 = −X,−Z)
    cam_el     = 1.55            # camera height above mid-plane (× lim)
    r_xz       = 3.96 * lim      # camera radius in the XZ plane
    az = np.radians(cam_az_deg)
    pl.camera_position = [
        (r_xz*np.cos(az), lim*cam_el, r_xz*np.sin(az)),  # camera position
        (0, 0, 0),                                        # focal point (target)
        (0, 1, 0),                                        # up vector (beam axis)
    ]

    if interactive:
        pl.show()   # opens live VTK window; rotate/zoom with mouse
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
    print(f"MM stack depth         : {t_MM*10:.2f} mm")
    print(f"PCB stack depth        : {t_PCB*10:.2f} mm")
    print(f"Trigger scint depth    : {t_scint*10:.2f} mm")
    print(f"LS stack depth         : {t_LS*10:.1f} mm  (3×{tLSCfrp*10:.0f}mm CFRP + 2×[{tLSInnerCfrp*1e3:.0f}µm CFRP + {tLSInnerAl*1e3:.0f}µm Al + {tLS*10:.0f}mm LAB])")
    print(f"Back scint depth       : {t_bsc*10:.1f} mm  (2cm plastic + 2×200µm tape)")
    print(f"Total stack depth      : {stack_depth*10:.1f} mm")
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
