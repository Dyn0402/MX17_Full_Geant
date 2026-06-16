#!/usr/bin/env python3
"""
compare_capsule_gas.py
Visual comparison of the He-3 gas region in the MX17 target capsule:
  OLD  — symmetric capsule used in the simulation until now
         (r=10 mm cylinder ±20 mm + two r=10 mm hemispheres, stops at ±30 mm)
  NEW  — full STEP-derived cavity: bore + lower hemisphere + conical fill
         channel up the neck to the Ø1.5 mm valve bore

Profiles are meridian (axis = beam = Y), drawn mirrored about the axis.
Run:  python scripts/compare_capsule_gas.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Profiles [mm] (z = beam axis, ro = radius) ───────────────────────────────
# OLD gas (DetectorConstruction.cc before the update)
Z_GAS_OLD = np.array([-30, -28, -26, -24, -22, -20,
                       20,  22,  24,  26,  28,  30.])
RO_GAS_OLD = np.array([0.001, 6.0, 8.0, 9.165, 9.798, 10.0,
                       10.0,  9.798, 9.165, 8.0, 6.0, 0.001])

# NEW gas (full STEP cavity)
Z_GAS_NEW = np.array([-29.5, -28, -26, -24, -22, -20, -15, -5, 5, 15, 20,
                      22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 44, 50.7])
RO_GAS_NEW = np.array([0.001, 6.0, 8.0, 9.165, 9.798, 10.0, 10.0, 10.0, 10.0,
                       10.0, 10.0, 9.798, 9.165, 8.0, 6.299, 4.842, 3.66,
                       2.711, 1.967, 1.41, 1.026, 0.75, 0.75])

# Al vessel outer wall (new 29-pt profile) — context only
Z_AL = np.array([-35, -34, -33, -31, -29, -27, -25, -23, -21, -20, -15, -5,
                 5, 15, 20, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 40, 45, 50, 51.])
RO_AL = np.array([0.0, 3.803, 5.287, 7.206, 8.48, 9.375, 9.994, 10.386, 10.6,
                  10.6, 10.6, 10.6, 10.6, 10.6, 10.6, 10.6, 10.386, 9.994,
                  9.375, 8.48, 7.206, 5.747, 4.708, 4.015, 3.621, 3.5, 3.5, 3.5, 3.5])


def mirror(z, ro):
    """Closed mirrored polygon (±ro, z) for a meridian profile."""
    return np.vstack([np.column_stack([ro, z]),
                      np.column_stack([-ro[::-1], z[::-1]])])


def draw_panel(ax, z_gas, ro_gas, title, ghost=None):
    # Al vessel wall (gray) behind the gas
    ax.add_patch(plt.Polygon(mirror(Z_AL, RO_AL), closed=True,
                             fc='#c8c8c8', ec='0.4', lw=0.8, zorder=2,
                             label='Al vessel wall'))
    # He-3 gas
    ax.add_patch(plt.Polygon(mirror(z_gas, ro_gas), closed=True,
                             fc='#4aa3df', ec='#1f6699', lw=1.0, zorder=3,
                             label='He-3 gas'))
    # optional dashed ghost of the other gas outline for reference
    if ghost is not None:
        gz, gro, glab = ghost
        ax.add_patch(plt.Polygon(mirror(gz, gro), closed=True,
                                 fc='none', ec='#d62728', lw=1.4, ls='--',
                                 zorder=4, label=glab))
    ax.axvline(0, color='0.6', lw=0.5, zorder=1)
    ax.set_aspect('equal')
    ax.set_xlim(-16, 16)
    ax.set_ylim(-38, 54)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel('radius  [mm]')
    ax.grid(True, lw=0.3, alpha=0.5)


fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 8), sharey=True)

draw_panel(ax0, Z_GAS_OLD, RO_GAS_OLD,
           'OLD (simulated)\nsymmetric capsule, gas stops at ±30 mm')
ax0.set_ylabel('Y  (beam axis)  [mm]')
ax0.legend(loc='lower right', fontsize=8, framealpha=0.9)

# New panel: show the old gas outline dashed so the added He-3 is obvious
draw_panel(ax1, Z_GAS_NEW, RO_GAS_NEW,
           'NEW (full STEP cavity)\ngas fills the neck fill-channel',
           ghost=(Z_GAS_OLD, RO_GAS_OLD, 'old gas extent'))
ax1.legend(loc='lower right', fontsize=8, framealpha=0.9)

# Annotate the added region on the new panel
ax1.annotate('added He-3\n(neck fill channel\nto valve bore)',
             xy=(2.5, 40), xytext=(8.5, 30), fontsize=8, color='#1f6699',
             ha='left', va='center',
             arrowprops=dict(arrowstyle='->', color='#1f6699', lw=1.0))

fig.suptitle('MX17 He-3 Target Capsule — gas region: OLD vs NEW',
             fontsize=13, fontweight='bold')
fig.tight_layout(rect=(0, 0, 1, 0.97))

_here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(_here, 'capsule_gas_old_vs_new.png')
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved: {out}")

# Quick numeric note on the He-3 gas volume change (solid of revolution)
def revolve_volume(z, ro):
    """Volume of a polycone (sum of conical frustа) [mm^3]."""
    V = 0.0
    for i in range(len(z) - 1):
        dz = z[i+1] - z[i]
        r1, r2 = ro[i], ro[i+1]
        V += np.pi * dz / 3.0 * (r1*r1 + r1*r2 + r2*r2)
    return V

v_old = revolve_volume(Z_GAS_OLD, RO_GAS_OLD)
v_new = revolve_volume(Z_GAS_NEW, RO_GAS_NEW)
print(f"He-3 gas volume  old: {v_old/1e3:6.3f} cm^3")
print(f"He-3 gas volume  new: {v_new/1e3:6.3f} cm^3  "
      f"(+{(v_new-v_old)/v_old*100:.1f} %)")
