#!/usr/bin/env python3
"""
plot_buildup.py — MX17 inside-out detector build-up (top-down slides).

Renders a 4-frame sequence showing the detector assembled from the INSIDE OUT,
for "building the detector" presentation slides:
  1. Micromegas (MM) only
  2. + SiPM trigger-scintillator walls
  3. + Liquid-scintillator (LS) layers
  4. + Back plastic-scintillator (PMT) bars     [full detector]

Two visual styles (choose with --style; default renders BOTH):
  clean    — flat blocks, coordinate axes, compass, big A/B/C/D labels
             (matches scripts/plot_mm_layout.py)
  detailed — every sub-layer shown (PCB, LS CFRP walls + 2 LAB layers,
             twin PMT bars with gap), richer legend
             (matches scripts/plot_geometry.py palette)

Coordinate convention (adopted 2026-06-30; see GEOMETRY_COORDINATE_CONVENTION.md):
  true top-down, +Z → right, +X → up, beam +Y ⊙ out of the page (right-handed).

Geometry (measured 2026-06-30; pinwheel shifts halved 2026-07-14):
  • Opposing mylar (window) faces: B↔D = 40.8 cm, C↔A = 40.9 cm.  Beam/target at
    the centre of the (roughly square) mylar box ⇒ faces at ±20.40 / ±20.45 cm.
  • Per-MM pinwheel (tangential) shift: D=1.55, B=1.575, A=1.635, C=1.73 cm.
  Layer thicknesses are imported from plot_geometry.py (SimConfig-derived).

Run:  python scripts/plot_buildup.py                 # both styles, 4 frames each
      python scripts/plot_buildup.py --style clean    # one style only
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon

# Layer depths / thicknesses come straight from plot_geometry.py (SimConfig),
# so the build-up stays in sync with the simulation stack.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_geometry as G   # noqa: E402

# ─────────────────────────────────────────────────────────────────────────────
# Geometry (measured; NEW spacing + halved pinwheel shift)
# ─────────────────────────────────────────────────────────────────────────────
SPAN_BD, SPAN_CA = 40.8, 40.9         # opposing mylar-face spans [cm]

MMS = [
    dict(letter='D', coord='+X', card='North', n=(+1.0, 0.0), dist=SPAN_BD/2, shift=1.55),
    dict(letter='B', coord='−X', card='South', n=(-1.0, 0.0), dist=SPAN_BD/2, shift=1.575),
    dict(letter='A', coord='+Z', card='East',  n=(0.0, +1.0), dist=SPAN_CA/2, shift=1.635),
    dict(letter='C', coord='−Z', card='West',  n=(0.0, -1.0), dist=SPAN_CA/2, shift=1.73),
]

def S(px, pz):
    """physics (X, Z) → screen (x = Z, y = X)."""
    return (pz, px)

# Face-relative layer-depth boundaries (from plot_geometry / SimConfig) [cm]
wMMf,  wMMb  = G.w_MM_f,  G.w_MM_b        # Micromegas stack (incl. drift)
wPCBf, wPCBb = G.w_PCB_f, G.w_PCB_b       # readout PCB
wSCf,  wSCb  = G.w_sc_f,  G.w_sc_b        # SiPM trigger scint wall
wLSf,  wLSb  = G.w_LS_f,  G.w_LS_b        # liquid-scintillator stack
wBSf,  wBSb  = G.w_bsc_f, G.w_bsc_b       # back plastic (PMT) scint bars
# LS internal sub-layers (CFRP walls + 2 LAB layers)
LS_SUBS = [
    (G._lsDivFront_f, G._lsDivFront_b, 'cfrp'),
    (G._wLS1_f,       G._wLS1_b,       'lab'),
    (G._lsDivMid_f,   G._lsDivMid_b,   'cfrp'),
    (G._wLS2_f,       G._wLS2_b,       'lab'),
    (G._lsDivRear_f,  G._lsDivRear_b,  'cfrp'),
]
HU = dict(mm=G.HW_U['mm'], scint=G.HW_U['scint'], ls=G.HW_U['ls'])
BS_HU, BS_OFF = G.bscTape_hu, G.bsc_u_offset

# He-3 capsule cross-section radii [cm]
HE3 = [(1.15, '#404040', 3), (1.06, '#b0b0b0', 4), (1.00, '#99d8f5', 5)]

# Palette
C_MM, C_PCB, C_SC = '#4a90d9', '#5cb85c', '#f0c040'
C_LS, C_CFRP, C_BS = '#d9534f', '#303030', '#e07820'
C_MYLAR = '#c0392b'

STAGE_NAME = {
    1: 'Micromegas (MM)',
    2: '+ SiPM trigger-scintillator walls',
    3: '+ Liquid-scintillator layers',
    4: '+ Back plastic-scintillator (PMT) bars   —   full detector',
}

# outermost radius (mylar face + full stack) → fixed frame extent for all frames
R_OUT = max(SPAN_BD, SPAN_CA) / 2.0 + wBSb
LIM = R_OUT + 6.0


# ─────────────────────────────────────────────────────────────────────────────
# Slab drawing
# ─────────────────────────────────────────────────────────────────────────────
def _slab(ax, arm, w_f, w_b, hu, color, u0=0.0, alpha=0.85, lw=0.6, z=2, ec='k'):
    """One layer slab for an arm, spanning depth [w_f, w_b] (from mylar face)
    and ±hu about the arm's tangential centre (shift + u0)."""
    n = np.array(arm['n'])
    t = np.array([-n[1], n[0]])
    u = arm['shift'] + u0
    d = arm['dist']
    pts = [n*(d+w_f) + t*(u+hu), n*(d+w_f) + t*(u-hu),
           n*(d+w_b) + t*(u-hu), n*(d+w_b) + t*(u+hu)]
    ax.add_patch(Polygon([S(*p) for p in pts], closed=True, facecolor=color,
                         edgecolor=ec, lw=lw, alpha=alpha, zorder=z))


def _mylar_front(ax, arm, hu):
    """Red MM entrance-window line on the front (target-facing) face."""
    n = np.array(arm['n'])
    t = np.array([-n[1], n[0]])
    u = arm['shift']
    d = arm['dist']
    p_p = n*d + t*(u+hu)
    p_m = n*d + t*(u-hu)
    ax.plot(*zip(S(*p_p), S(*p_m)), color=C_MYLAR, lw=2.5, zorder=6)


def _draw_arm_layers(ax, arm, stage, style):
    """Draw all layers up to `stage` for one arm, in the given style."""
    detailed = (style == 'detailed')

    # Stage 1 — Micromegas (+ PCB readout in detailed style)
    _slab(ax, arm, wMMf, wMMb, HU['mm'], C_MM, alpha=0.88, z=4)
    if detailed:
        _slab(ax, arm, wPCBf, wPCBb, HU['mm'], C_PCB, alpha=0.85, z=4)
    _mylar_front(ax, arm, HU['mm'])
    if stage < 2:
        return

    # Stage 2 — SiPM trigger-scintillator wall
    _slab(ax, arm, wSCf, wSCb, HU['scint'], C_SC, alpha=0.85, z=3)
    if stage < 3:
        return

    # Stage 3 — Liquid scintillator
    if detailed:
        for lf, lb, kind in LS_SUBS:
            _slab(ax, arm, lf, lb, HU['ls'],
                  C_CFRP if kind == 'cfrp' else C_LS,
                  alpha=0.85 if kind == 'cfrp' else 0.75, z=3)
    else:
        _slab(ax, arm, wLSf, wLSb, HU['ls'], C_LS, alpha=0.8, z=3)
    if stage < 4:
        return

    # Stage 4 — Back plastic (PMT) scint bars: two per arm with a gap
    for sgn in (-1.0, +1.0):
        _slab(ax, arm, wBSf, wBSb, BS_HU, C_BS, u0=sgn*BS_OFF, alpha=0.85, z=3)


# ─────────────────────────────────────────────────────────────────────────────
# Shared decorations (identical across every frame → slides overlay cleanly)
# ─────────────────────────────────────────────────────────────────────────────
def _decorate(ax, stage, style):
    # He-3 target capsule
    for r, fc, zo in HE3:
        ax.add_patch(plt.Circle((0, 0), r, color=fc, zorder=zo))

    # Arm labels — fixed radius (just beyond the full stack) so they never jump
    for arm in MMS:
        n = np.array(arm['n'])
        t = np.array([-n[1], n[0]])
        pos = n*(arm['dist'] + wBSb + 3.2) + t*arm['shift']
        sx, sy = S(*pos)
        ax.text(sx, sy, f"MM {arm['letter']}\n{arm['coord']} · {arm['card']}",
                ha='center', va='center', fontsize=11.5, fontweight='bold',
                linespacing=1.35, zorder=8,
                bbox=dict(boxstyle='round,pad=0.32', facecolor='white',
                          edgecolor=C_MM, alpha=0.95))

    # Coordinate axes through origin: +Z right, +X up
    L = 11.0
    ax.annotate('', xy=S(*(0,  L)), xytext=(0, 0),
                arrowprops=dict(arrowstyle='-|>', color='k', lw=2.0), zorder=7)
    ax.annotate('', xy=S(*(L, 0)), xytext=(0, 0),
                arrowprops=dict(arrowstyle='-|>', color='k', lw=2.0), zorder=7)
    ax.plot([0, S(*(0, -L))[0]], [0, S(*(0, -L))[1]], color='0.6', lw=1.0, ls='--', zorder=1)
    ax.plot([0, S(*(-L, 0))[0]], [0, S(*(-L, 0))[1]], color='0.6', lw=1.0, ls='--', zorder=1)
    ax.text(*S(*(0, L + 1.3)), '+Z', ha='left',   va='center', fontsize=12, fontweight='bold', zorder=7)
    ax.text(*S(*(L + 1.3, 0)), '+X', ha='center', va='bottom', fontsize=12, fontweight='bold', zorder=7)

    # Beam: +Y out of the page (⊙)
    ax.scatter([0], [0], s=420, facecolors='none', edgecolors='firebrick', linewidths=2.0, zorder=9)
    ax.scatter([0], [0], s=40, color='firebrick', zorder=9)
    ax.text(-1.8, -1.8, 'beam +Y\n(out of page)', color='firebrick',
            fontsize=8.5, ha='right', va='top', zorder=9)

    # Compass rose (top-right)
    cx, cy, cr = LIM - 6.5, LIM - 6.5, 4.0
    for ang, lab in [(90, 'N'), (0, 'E'), (270, 'S'), (180, 'W')]:
        dx, dy = cr*np.cos(np.radians(ang)), cr*np.sin(np.radians(ang))
        ax.annotate('', xy=(cx+dx, cy+dy), xytext=(cx, cy),
                    arrowprops=dict(arrowstyle='-|>', color='#2c6e49', lw=1.4), zorder=8)
        ax.text(cx + 1.45*dx, cy + 1.45*dy, lab, ha='center', va='center',
                color='#2c6e49', fontsize=10.5, fontweight='bold', zorder=8)
    ax.text(cx, cy - cr - 3.6, '+X=N · +Z=E', ha='center', va='center',
            color='#2c6e49', fontsize=7.5, zorder=8)

    # Legend (cumulative up to this stage)
    handles = [mpatches.Patch(color=C_MM, alpha=0.88,
                              label='Micromegas' + (' (drift)' if style == 'detailed' else ''))]
    if style == 'detailed':
        handles.append(mpatches.Patch(color=C_PCB, alpha=0.85, label='Readout PCB'))
    handles.append(mpatches.Patch(color=C_MYLAR, label='MM mylar window (faces target)'))
    if stage >= 2:
        handles.append(mpatches.Patch(color=C_SC, alpha=0.85,
                                      label='SiPM trigger scint  (48×48 cm)'))
    if stage >= 3:
        handles.append(mpatches.Patch(color=C_LS, alpha=0.78,
                                      label='Liquid scint LAB  (2×2 cm, 45×45 cm)'))
        if style == 'detailed':
            handles.append(mpatches.Patch(color=C_CFRP, alpha=0.85, label='LS CFRP walls + liners'))
    if stage >= 4:
        handles.append(mpatches.Patch(color=C_BS, alpha=0.85,
                                      label='Back plastic scint / PMT bars  (2×[20×30] cm)'))
    handles.append(mpatches.Patch(color='#99d8f5', label='He-3 target gas bore (Ø20 mm)'))
    ax.legend(handles=handles, loc='lower left', fontsize=8, framealpha=0.92)


def make_frame(stage, style):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')
    for arm in MMS:
        _draw_arm_layers(ax, arm, stage, style)
    _decorate(ax, stage, style)

    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(-LIM, LIM)
    ax.set_xlabel('Z  [cm]  (East →)', fontsize=11)
    ax.set_ylabel('X  [cm]  (North ↑)', fontsize=11)
    ax.set_title(f'MX17 Detector Build-Up   ·   Step {stage} / 4\n'
                 f'{STAGE_NAME[stage]}\n'
                 f'+Z right · +X up · beam +Y ⊙ out of page   ({style} style)',
                 fontsize=11)
    ax.axhline(0, color='0.88', lw=0.5, zorder=0)
    ax.axvline(0, color='0.88', lw=0.5, zorder=0)
    ax.grid(True, lw=0.3, alpha=0.4)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    styles = ['clean', 'detailed']
    if '--style' in sys.argv:
        styles = [sys.argv[sys.argv.index('--style') + 1]]

    tags = {1: '1_mm', 2: '2_sipm', 3: '3_ls', 4: '4_full'}
    print(f"Outer stack radius : {R_OUT:.1f} cm   (frame ±{LIM:.0f} cm)")
    for style in styles:
        for stage in (1, 2, 3, 4):
            fig = make_frame(stage, style)
            base = os.path.join(here, f'mx17_buildup_{style}_{tags[stage]}')
            fig.savefig(base + '.png', dpi=150, bbox_inches='tight')
            fig.savefig(base + '.pdf', bbox_inches='tight')
            plt.close(fig)
            print(f"Saved: {os.path.basename(base)}.png/.pdf")
