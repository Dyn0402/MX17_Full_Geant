#!/usr/bin/env python3
"""
plot_mm_layout.py
Top-down experimental layout of the four MX17 Micromegas (MM only — no
scintillators yet).

ADOPTED EXPERIMENT CONVENTION (true top-down, right-handed)
  • View      : looking straight DOWN from above (bird's-eye).
  • Screen    : +Z → right, +X → up   (i.e. screen_x = Z, screen_y = X)
  • Beam +Y   : OUT of the page toward the viewer (⊙).  Right-handed:
                ŷ = ẑ × x̂ = (right) × (up) = out of page, so +Y = beam,
                consistent with the Geant sim (SetParticleMomentumDirection
                (0,1,0)).  The map sim(X→up, Y→out, Z→right) is a cyclic
                axis permutation = a pure rotation (no mirror).

  MM labels & cardinals (user-defined):
      +X = D = North (up)        −X = B = South (down)
      +Z = A = East  (right)     −Z = C = West  (left)

GEOMETRY UPDATE 2026-06-30 (measured at the experiment, refined; NOT yet
propagated to the Geant sim / SimConfig — to be fixed everywhere later):
  • Opposing mylar (window) face spacing:  B↔D = 40.8 cm,  C↔A = 40.9 cm.
  • Beam/target placed at the CENTRE of the (roughly square) mylar-face box
    ⇒ mylar faces at ±20.40 cm (X pair) and ±20.45 cm (Z pair).
  • Per-MM pinwheel shift (tangential): D=1.55, A=1.635, B=1.575, C=1.73 cm
    (halved 2026-07-14 — the earlier 3.10–3.46 cm offsets were 2× too large;
    no longer a uniform 30 mm).
MM active width 38 cm; stack depth ≈ 3.04 cm (incl. 30 mm drift gap).
Run:  python scripts/plot_mm_layout.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon

# ─────────────────────────────────────────────────────────────────────────────
# Geometry (cm) — from SimConfig / plot_geometry.py
# ─────────────────────────────────────────────────────────────────────────────
# ── Opposing mylar-face spacing (measured 2026-06-30, refined) ───────────────
# Distances between the inner mylar (entrance-window) faces of opposing MMs.
# The set-up is positioned so the target/beam sit at the CENTRE of the (roughly
# square) box formed by these four mylar faces ⇒ each face is half its pair-span
# from the beam axis.
SPAN_BD = 40.8                    # B(−X) ↔ D(+X) mylar-face distance [cm]  → ±20.40
SPAN_CA = 40.9                    # C(−Z) ↔ A(+Z) mylar-face distance [cm]  → ±20.45
MM_U_HALF = 38.0 / 2.0           # MM in-plane half-width (active u = 38 cm)

# ── Per-MM pinwheel (circular) tangential shift ──────────────────────────────
# Each MM is shifted tangentially (⟂ its outward normal) by a measured amount,
# in a right-handed pinwheel: shift = +t̂ · |shift|, with t̂ = (−n_z, n_x).
# Directions (screen): B→left(−Z), C→up(+X), D→right(+Z), A→down(−X).
# Magnitudes are now individual (no longer a uniform 30 mm).

# MM stack depth (mylar + windows + kapton/Cu cathode + 30 mm drift + mesh +
# amplification + resistive paste), same sum as plot_geometry.py:
_um, _mm = 1e-4, 0.1
MM_DEPTH = (40*_um + 0.1*_um + 50*_um + 9*_um + 30*_mm
            + 30*_um + 150*_um + 100*_um)          # ≈ 3.04 cm

# He-3 target capsule (top-down cross-section)
HE3_GAS_R  = 1.00         # gas bore radius [cm]
HE3_AL_R   = 1.06         # Al outer radius [cm]
HE3_CFRP_R = 1.15         # CFRP outer radius [cm]

# ─────────────────────────────────────────────────────────────────────────────
# Detector definitions.  n = outward normal in physics (X, Z).
# ─────────────────────────────────────────────────────────────────────────────
# dist = mylar-face distance from beam axis (= half the opposing pair-span,
#        so the target sits at the centre of the mylar box).
# shift = measured tangential pinwheel offset [cm].
MMS = [
    dict(letter='D', coord='+X', card='North', n=(+1.0, 0.0), dist=SPAN_BD/2, shift=1.55),
    dict(letter='B', coord='−X', card='South', n=(-1.0, 0.0), dist=SPAN_BD/2, shift=1.575),
    dict(letter='A', coord='+Z', card='East',  n=(0.0, +1.0), dist=SPAN_CA/2, shift=1.635),
    dict(letter='C', coord='−Z', card='West',  n=(0.0, -1.0), dist=SPAN_CA/2, shift=1.73),
]

# screen mapping: physics (X, Z) → screen (x=Z, y=X)
def S(px, pz):
    return (pz, px)


def plot_mm_layout():
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')

    # ── He-3 target capsule at origin (gas / Al / CFRP) ──────────────────────
    for r, fc, zo in [(HE3_CFRP_R, '#404040', 3),
                      (HE3_AL_R,   '#b0b0b0', 4),
                      (HE3_GAS_R,  '#99d8f5', 5)]:
        ax.add_patch(plt.Circle((0, 0), r, color=fc, zorder=zo))
    ax.annotate('He-3 target',
                xy=(0, 0), xytext=(HE3_CFRP_R + 1.0, HE3_CFRP_R + 1.0),
                fontsize=8, zorder=6,
                arrowprops=dict(arrowstyle='-', lw=0.6, color='0.4'))

    # ── Four Micromegas ──────────────────────────────────────────────────────
    for d in MMS:
        n = np.array(d['n'])               # outward normal (X, Z)
        t = np.array([-n[1], n[0]])        # in-plane tangent ⟂ n
        dist = d['dist']                   # mylar-face distance from axis
        shift = t * d['shift']             # pinwheel tangential offset
        c_f = n * dist + shift             # mylar-face centre (X, Z), shifted

        # corners (X, Z): front ±t, then back face (+depth along n)
        cf_p = c_f + t * MM_U_HALF
        cf_m = c_f - t * MM_U_HALF
        cb_m = c_f + n * MM_DEPTH - t * MM_U_HALF
        cb_p = c_f + n * MM_DEPTH + t * MM_U_HALF
        poly = [S(*cf_p), S(*cf_m), S(*cb_m), S(*cb_p)]
        ax.add_patch(Polygon(poly, closed=True, facecolor='#4a90d9',
                             edgecolor='k', lw=1.0, alpha=0.85, zorder=2))

        # entrance mylar window (front face, toward target) highlighted in red
        ax.plot(*zip(S(*cf_p), S(*cf_m)), color='#c0392b', lw=2.5, zorder=3)

        # pinwheel shift indicator: arrow from un-shifted to shifted face centre
        usc, ssc = S(*(n * dist)), S(*c_f)
        ax.annotate('', xy=ssc, xytext=usc,
                    arrowprops=dict(arrowstyle='-|>', color='#7b4fa3', lw=1.6),
                    zorder=7)
        # magnitude label, nudged inward (toward centre) into open space
        mlx, mly = S(*(n * dist + shift / 2.0 - n * 1.7))
        _slab = f"{d['shift']:.3f}".rstrip('0').rstrip('.')
        ax.text(mlx, mly, f"{_slab} cm", color='#7b4fa3', fontsize=7.5,
                ha='center', va='center', zorder=7)

        # single boxed label, nudged outward so it sits clear of the slab
        cen = c_f + n * (MM_DEPTH / 2.0)
        off = n * 6.5
        sx, sy = S(*(cen + off))
        ax.text(sx, sy, f"MM {d['letter']}\n{d['coord']} · {d['card']}",
                ha='center', va='center', fontsize=12, fontweight='bold',
                linespacing=1.4, zorder=6,
                bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                          edgecolor='#4a90d9', alpha=0.95))

    # ── Opposing mylar-face dimensions (B↔D and C↔A) ─────────────────────────
    hx, hz = SPAN_BD / 2.0, SPAN_CA / 2.0      # face half-spans (X, Z)
    # B↔D : vertical (along X), drawn at screen x = Zdim
    Zdim = -15.0
    ax.annotate('', xy=(Zdim, hx), xytext=(Zdim, -hx),
                arrowprops=dict(arrowstyle='<|-|>', color='#7b4fa3', lw=1.4), zorder=4)
    ax.text(Zdim - 0.9, 0, f'B↔D  {SPAN_BD:.1f} cm', color='#7b4fa3', fontsize=8.5,
            rotation=90, ha='right', va='center', fontweight='bold', zorder=4)
    # C↔A : horizontal (along Z), drawn at screen y = Xdim
    Xdim = -15.0
    ax.annotate('', xy=(hz, Xdim), xytext=(-hz, Xdim),
                arrowprops=dict(arrowstyle='<|-|>', color='#7b4fa3', lw=1.4), zorder=4)
    ax.text(0, Xdim - 0.9, f'C↔A  {SPAN_CA:.1f} cm', color='#7b4fa3', fontsize=8.5,
            ha='center', va='top', fontweight='bold', zorder=4)

    # ── Coordinate axes through origin: +Z right, +X up ──────────────────────
    L = 13.0
    ax.annotate('', xy=S(*(0,  L)), xytext=(0, 0),   # +Z (screen right)
                arrowprops=dict(arrowstyle='-|>', color='k', lw=2.0))
    ax.annotate('', xy=S(*( L, 0)), xytext=(0, 0),   # +X (screen up)
                arrowprops=dict(arrowstyle='-|>', color='k', lw=2.0))
    # faint negative axes
    ax.plot([0, S(*(0, -L))[0]], [0, S(*(0, -L))[1]], color='0.6', lw=1.0, ls='--')
    ax.plot([0, S(*(-L, 0))[0]], [0, S(*(-L, 0))[1]], color='0.6', lw=1.0, ls='--')
    ax.text(*S(*(0, L + 1.2)), '+Z', ha='left',   va='center', fontsize=12, fontweight='bold')
    ax.text(*S(*(L + 1.2, 0)), '+X', ha='center', va='bottom', fontsize=12, fontweight='bold')

    # ── Beam: +Y out of the page (true top-down) ─────────────────────────────
    ax.scatter([0], [0], s=420, facecolors='none', edgecolors='firebrick',
               linewidths=2.0, zorder=7)
    ax.scatter([0], [0], s=40, color='firebrick', zorder=7)   # the dot of ⊙
    ax.text(-2.0, -2.0, 'beam +Y\n(out of page)\n@ mylar-box centre',
            color='firebrick', fontsize=9, ha='right', va='top', zorder=7)

    # ── Cardinal compass rose (inset, top-right corner) ──────────────────────
    cxx, cyy, cr = 33.0, 33.0, 4.5     # compass centre & radius (screen)
    for ang, lab in [(90, 'N'), (0, 'E'), (270, 'S'), (180, 'W')]:
        dx, dy = cr*np.cos(np.radians(ang)), cr*np.sin(np.radians(ang))
        ax.annotate('', xy=(cxx+dx, cyy+dy), xytext=(cxx, cyy),
                    arrowprops=dict(arrowstyle='-|>', color='#2c6e49', lw=1.5),
                    zorder=8)
        ax.text(cxx + 1.45*dx, cyy + 1.45*dy, lab, ha='center', va='center',
                color='#2c6e49', fontsize=11, fontweight='bold', zorder=8)
    ax.text(cxx, cyy - cr - 2.8, '+X=N  ·  +Z=E', ha='center', va='center',
            color='#2c6e49', fontsize=8, zorder=8)

    # ── Legend ───────────────────────────────────────────────────────────────
    handles = [
        mpatches.Patch(color='#4a90d9', alpha=0.85,
                       label=f'Micromegas  ({2*MM_U_HALF:.0f} cm wide, '
                             f'{MM_DEPTH:.1f} cm deep)'),
        mpatches.Patch(color='#c0392b', label='MM entrance mylar window (faces target)'),
        mpatches.Patch(color='#99d8f5', label='He-3 target gas bore (Ø20 mm)'),
    ]
    ax.legend(handles=handles, loc='lower left', fontsize=8, framealpha=0.9)

    # ── Frame ────────────────────────────────────────────────────────────────
    lim = max(SPAN_BD, SPAN_CA) / 2.0 + MM_DEPTH + 14
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel('Z  [cm]  (East →)', fontsize=11)
    ax.set_ylabel('X  [cm]  (North ↑)', fontsize=11)
    ax.set_title('MX17 Micromegas — Top-Down Layout (beam @ mylar-box centre)\n'
                 '+Z right · +X up · beam +Y ⊙ out of page (right-handed)\n'
                 'mylar faces B↔D 40.8 cm · C↔A 40.9 cm · pinwheel shift 1.55–1.73 cm',
                 fontsize=10.5)
    ax.axhline(0, color='0.85', lw=0.5, zorder=0)
    ax.axvline(0, color='0.85', lw=0.5, zorder=0)
    ax.grid(True, lw=0.3, alpha=0.5)
    fig.tight_layout()
    return fig


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, 'mx17_mm_layout_topdown')
    print(f"B↔D mylar span         : {SPAN_BD:.1f} cm  → faces at ±{SPAN_BD/2:.2f} cm")
    print(f"C↔A mylar span         : {SPAN_CA:.1f} cm  → faces at ±{SPAN_CA/2:.2f} cm")
    print(f"Beam/target            : at centre of mylar box")
    print("Pinwheel shifts        : " +
          ", ".join(f"{d['letter']}={d['shift']:.2f}cm" for d in MMS))
    print(f"MM width (active u)    : {2*MM_U_HALF:.1f} cm")
    print(f"MM stack depth         : {MM_DEPTH:.2f} cm")
    fig = plot_mm_layout()
    fig.savefig(out + '.png', dpi=150, bbox_inches='tight')
    fig.savefig(out + '.pdf', bbox_inches='tight')
    print(f"Saved: {out}.png/.pdf")
