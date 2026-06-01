#!/usr/bin/env python3
"""
plot_geometry.py
Visualise the MX17 4-arm detector geometry.
  Figure 1 — top-down 2D (XZ plane, beam/Y axis into the page)  [matplotlib]
  Figure 2 — 3D isometric view                                   [pyvista/VTK]

No Geant4 dependency; dimensions are taken directly from SimConfig defaults.
Run from anywhere:  python scripts/plot_geometry.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import pyvista as pv

pv.OFF_SCREEN = True   # headless rendering — no display required


# ─────────────────────────────────────────────────────────────────────────────
# SimConfig defaults  (keep in sync with include/SimConfig.hh)
# ─────────────────────────────────────────────────────────────────────────────
CFG = dict(
    mm_distance_cm        = 22.0,
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
    gap_ls_to_back_mm     = 10.0,
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


# ─────────────────────────────────────────────────────────────────────────────
# Arm geometry
# ─────────────────────────────────────────────────────────────────────────────
dist = CFG['mm_distance_cm']

ARM_DEF = [
    dict(id=0, label='+X', ff=np.array([ dist, 0,    0   ]),
         u_hat=np.array([0,0,-1]), w_hat=np.array([ 1,0,0])),
    dict(id=1, label='−X', ff=np.array([-dist, 0,    0   ]),
         u_hat=np.array([0,0, 1]), w_hat=np.array([-1,0,0])),
    dict(id=2, label='+Z', ff=np.array([   0,  0,  dist  ]),
         u_hat=np.array([1,0, 0]), w_hat=np.array([0,0, 1])),
    dict(id=3, label='−Z', ff=np.array([   0,  0, -dist  ]),
         u_hat=np.array([-1,0,0]), w_hat=np.array([0,0,-1])),
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
    ('MM',       w_MM_f,  w_MM_b,  'mm',    'mm',    '#4a90d9'),
    ('PCB',      w_PCB_f, w_PCB_b, 'mm',    'mm',    '#5cb85c'),
    ('PlScint',  w_sc_f,  w_sc_b,  'scint', 'scint', '#f0c040'),
    ('LiqScint', w_LS_f,  w_LS_b,  'ls',    'ls',    '#d9534f'),
]

he3_r       = 1.5   # radius = 1.5 cm  → diameter = 3 cm
al_wall     = 0.05
cfrp_wall   = 0.09
he3_total_r = he3_r + al_wall + cfrp_wall
he3_half_y  = 4.0 + al_wall + cfrp_wall  # half-length = 4 cm → length = 8 cm


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

    # He-3 capsule cross-section — largest drawn first (behind), gas on top
    for r, fc, zo in [
        (he3_total_r,     '#404040', 3),
        (he3_r + al_wall, '#b0b0b0', 4),
        (he3_r,           '#99d8f5', 5),
    ]:
        ax.add_patch(plt.Circle((0, 0), r, color=fc, zorder=zo))

    # Arm detector slabs projected onto XZ
    for arm in ARM_DEF:
        for lname, w_f, w_b, ukey, vkey, col in LAYERS_2D:
            cen, hs = arm_box_world(arm, w_f, w_b, HW_U[ukey], 0)
            x0, z0 = cen[0] - hs[0], cen[2] - hs[2]
            ax.add_patch(Rectangle((x0, z0), 2*hs[0], 2*hs[2],
                                    linewidth=0.5, edgecolor='k',
                                    facecolor=col, alpha=0.75, zorder=2))
        ff  = arm['ff']
        mid = ff + arm['w_hat'] * (stack_depth / 2)
        loff = ff / np.linalg.norm(ff) * 2.5
        ax.text(mid[0]+loff[0], mid[2]+loff[2],
                f"Arm {arm['id']}\n({arm['label']})",
                ha='center', va='center', fontsize=9, fontweight='bold', zorder=5)

    # Back plastic scints — two individual bars per arm with gap visible
    for arm in ARM_DEF:
        for u_sign in [-1.0, +1.0]:
            cen, hs = arm_box_world(arm, w_bsc_f, w_bsc_b,
                                    bscTape_hu, bscTape_hv,
                                    u_offset=u_sign * bsc_u_offset)
            x0, z0 = cen[0] - hs[0], cen[2] - hs[2]
            ax.add_patch(Rectangle((x0, z0), 2*hs[0], 2*hs[2],
                                    linewidth=0.5, edgecolor='k',
                                    facecolor='#e07820', alpha=0.80, zorder=2))

    # Beam indicator (no dot — just the label)
    ax.annotate('beam ⊙\n(+Y)', xy=(0, 0), xytext=(3, 3), fontsize=8, zorder=6)

    legend_patches = [
        mpatches.Patch(color='#99d8f5', label=f'He-3 gas  (r = {he3_r} cm)'),
        mpatches.Patch(color='#404040', label='He-3 capsule walls (Al 0.5 mm + CFRP 0.9 mm)'),
        mpatches.Patch(color='#4a90d9', alpha=0.75, label=f'MM drift gas  ({tDrift*10:.0f} mm drift)'),
        mpatches.Patch(color='#5cb85c', alpha=0.75, label=f'PCB stack  ({t_PCB*10:.1f} mm)'),
        mpatches.Patch(color='#f0c040', alpha=0.75,
                       label=f'Trigger plastic scint  ({tPlScint*10:.0f} mm, 48×48 cm)'),
        mpatches.Patch(color='#d9534f', alpha=0.75,
                       label=f'Liq. scint. stack  (2×{tLS*10:.0f} mm LAB, {CFG["ls_size_u_cm"]:.0f}×{CFG["ls_size_v_cm"]:.0f} cm)'),
        mpatches.Patch(color='#e07820', alpha=0.75,
                       label=f'Back plastic scints  ({int(bsc_u)} mm, {int(bsc_v)}×(2×{int(bsc_th*10)}) cm)'),
    ]
    ax.legend(handles=legend_patches, loc='upper left', fontsize=8, framealpha=0.6)

    lim = dist + stack_depth + 4
    ax.set_xlim(-lim, lim);  ax.set_ylim(-lim, lim)
    ax.set_xlabel('X  [cm]', fontsize=12)
    ax.set_ylabel('Z  [cm]', fontsize=12)
    ax.set_title('MX17 Detector Geometry — Top-Down View (XZ plane)\n'
                 'Beam along +Y (into page)', fontsize=12)
    ax.axhline(0, color='0.7', lw=0.5, zorder=1)
    ax.axvline(0, color='0.7', lw=0.5, zorder=1)
    ax.grid(True, lw=0.3, alpha=0.5)

    _annotate_dim(ax, 0,    dist,
                  y=-(lim - 1.5), label=f'{dist:.0f} cm')
    _annotate_dim(ax, dist, dist + stack_depth,
                  y=-(lim - 1.5), label=f'{stack_depth*10:.0f} mm')

    fig.tight_layout()
    return fig


def _annotate_dim(ax, x1, x2, y, label, color='0.4'):
    ax.annotate('', xy=(x2, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle='<->', color=color, lw=1.2))
    ax.text((x1+x2)/2, y - 0.8, label,
            ha='center', va='top', fontsize=7, color=color)


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
        ('LS',  w_LS_f,  w_LS_b,  'ls',    (0.85, 0.33, 0.31), 0.70,
         f'Liq. scint. stack  (2×{tLS*10:.0f} mm LAB, {CFG["ls_size_u_cm"]:.0f}×{CFG["ls_size_v_cm"]:.0f} cm)'),
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

    # He-3 capsule: nested solid cylinders along Y (outer first so inner shows through)
    caps = [
        (he3_total_r,     (0.16, 0.16, 0.16), 0.70, 'CFRP capsule wall'),
        (he3_r + al_wall, (0.67, 0.67, 0.67), 0.70, 'Al capsule wall'),
        (he3_r,           (0.60, 0.85, 0.96), 0.85, f'He-3 gas  (r = {he3_r} cm)'),
    ]
    for r, col, alpha, llabel in caps:
        cyl = pv.Cylinder(center=(0, 0, 0), direction=(0, 1, 0),
                          radius=r, height=2*he3_half_y, resolution=80,
                          capping=True)
        pl.add_mesh(cyl, color=col, opacity=alpha,
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

    # Camera: from well outside the −X,+Y,−Z corner looking at origin; Y is "up"
    lim = dist + stack_depth + 5
    pl.camera_position = [
        (-lim * 2.8,  lim * 1.2, -lim * 2.8),  # camera position
        (0, 0, 0),                                # focal point
        (0, 1, 0),                                # up vector (beam axis)
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
    interactive = True
    _here = os.path.dirname(os.path.abspath(__file__))

    print(f"MM front face distance : {dist:.1f} cm")
    print(f"MM stack depth         : {t_MM*10:.2f} mm")
    print(f"PCB stack depth        : {t_PCB*10:.2f} mm")
    print(f"Trigger scint depth    : {t_scint*10:.2f} mm")
    print(f"LS stack depth         : {t_LS*10:.1f} mm  (3×{tLSCfrp*10:.0f}mm CFRP + 2×[{tLSInnerCfrp*1e3:.0f}µm CFRP + {tLSInnerAl*1e3:.0f}µm Al + {tLS*10:.0f}mm LAB])")
    print(f"Back scint depth       : {t_bsc*10:.1f} mm  (2cm plastic + 2×200µm tape)")
    print(f"Total stack depth      : {stack_depth*10:.1f} mm")
    print(f"Arm outer edge         : {dist + stack_depth:.1f} cm from origin")

    out1 = os.path.join(_here, 'mx17_geometry_topdown.png')
    out2 = os.path.join(_here, 'mx17_geometry_3d.png')

    fig1 = plot_2d_topdown()
    fig1.savefig(out1, dpi=150, bbox_inches='tight')
    print(f"Saved: {out1}")

    plot_3d_pyvista(out_path=out2, interactive=interactive)
    plt.show()
