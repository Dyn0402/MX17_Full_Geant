#!/usr/bin/env python3
"""
Plot EAR2 neutron beam spatial (radial) profile vs energy from the FLUKA transport simulation.
Data: lamda2DvsEn_EAR2.root  (from Marta Sabaté-Gilarte's AFS area)
  Lambda2D: TH2D, axes = [neutron energy (eV), radial distance from beam axis (cm)]
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import uproot
from pathlib import Path
from scipy.ndimage import uniform_filter1d

DATA_DIR = Path('/home/dylan/CLionProjects/MX17_Full_Geant/data')
OUT_DIR  = Path('/home/dylan/CLionProjects/MX17_Full_Geant/data')

PROFILE_FILE = DATA_DIR / 'lamda2DvsEn_EAR2.root'

# X17 ROI
E_X17_LOW_EV  = 10.0e6
E_X17_HIGH_EV = 20.0e6

# Radial rebin factor (original bin = 0.001 cm; rebin to 0.05 cm = 50x)
REBIN_R = 50

# Energy slices to show in radial profile plots [eV]
E_SLICES = {
    '1 meV':  1e-3,
    '1 eV':   1.0,
    '1 keV':  1e3,
    '100 keV': 1e5,
    '1 MeV':  1e6,
    '10 MeV': 1e7,
    '20 MeV': 2e7,
}

# Target radii to evaluate (for He-3 target context)
R_TARGETS_CM = [0.5, 1.0, 1.5, 2.0, 2.5]


def load_2d():
    with uproot.open(PROFILE_FILE) as f:
        h = f['Lambda2D']
        vals = h.values().copy()        # shape: (nE, nR)
        ex   = h.axis(0).edges().copy() # energy edges
        ey   = h.axis(1).edges().copy() # radial edges
    cx = 0.5 * (ex[:-1] + ex[1:])
    cy = 0.5 * (ey[:-1] + ey[1:])
    return vals, cx, ex, cy, ey


def rebin_r(vals, cy, ey, factor):
    nE, nR = vals.shape
    nR_new = nR // factor
    v_rb = vals[:, :nR_new*factor].reshape(nE, nR_new, factor).sum(axis=2)
    ey_rb = ey[:nR_new*factor+1:factor]
    cy_rb = 0.5 * (ey_rb[:-1] + ey_rb[1:])
    return v_rb, cy_rb, ey_rb


def nearest_energy_idx(cx, e_target):
    return int(np.argmin(np.abs(cx - e_target)))


def main():
    vals, cx, ex, cy_fine, ey_fine = load_2d()

    # Rebinned radial axis
    vals_rb, cy, ey = rebin_r(vals, cy_fine, ey_fine, REBIN_R)
    dr = ey[1] - ey[0]

    plot_2d_map(vals_rb, cx, ex, cy, ey)
    plot_radial_slices(vals_rb, cx, cy, dr)
    plot_beam_size_vs_energy(vals_rb, cx, cy)
    plot_target_fraction_vs_energy(vals_rb, cx, cy, ey)
    plt.show()


def plot_2d_map(vals_rb, cx, ex, cy, ey):
    # Normalize each energy slice to its max so the profile shape is visible
    # independent of the absolute flux
    norm_map = vals_rb.copy()
    col_max = norm_map.max(axis=1, keepdims=True)
    col_max[col_max == 0] = 1
    norm_map /= col_max

    fig, ax = plt.subplots(figsize=(11, 6))
    pcm = ax.pcolormesh(ex, ey, norm_map.T,
                        norm=mcolors.PowerNorm(gamma=0.4, vmin=0, vmax=1),
                        cmap='inferno', shading='flat', rasterized=True)
    cb = fig.colorbar(pcm, ax=ax, label='Relative intensity (normalised per energy slice)')
    ax.set_xscale('log')
    ax.set_xlabel('Neutron energy (eV)')
    ax.set_ylabel('Radial distance from beam axis (cm)')
    ax.set_ylim(0, 3.0)
    ax.set_xlim(ex[0], ex[-1])
    ax.axhline(2.5, color='cyan', linewidth=1.0, linestyle='--', label='Collimator r = 2.5 cm')
    ax.axvspan(E_X17_LOW_EV, E_X17_HIGH_EV, alpha=0.25, color='lime',
               label='X17 ROI 10–20 MeV')
    ax.legend(fontsize=9)
    ax.set_title('EAR2 Neutron Beam — Radial Profile vs Energy (FLUKA transport)')
    _add_energy_ticks(ax)
    fig.tight_layout()
    _save(fig, 'ear2_spatial_2d_map.png')


def plot_radial_slices(vals_rb, cx, cy, dr):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle('EAR2 Radial Beam Profile at Selected Energies', fontsize=13)

    colors = plt.cm.plasma(np.linspace(0.05, 0.95, len(E_SLICES)))

    for ax_idx, (ax, norm_by) in enumerate(zip(axes, ['peak', 'area'])):
        for (label, e_val), color in zip(E_SLICES.items(), colors):
            ie = nearest_energy_idx(cx, e_val)
            profile = vals_rb[ie, :].copy()

            # Smooth lightly (3-bin running average) for display
            profile_sm = uniform_filter1d(profile, size=3)

            if norm_by == 'peak':
                peak = profile_sm.max()
                if peak > 0:
                    profile_sm = profile_sm / peak
                ylabel = 'Intensity (normalised to peak)'
            else:
                # normalise to unit integral over r (probability density in r)
                total = profile_sm.sum() * dr
                if total > 0:
                    profile_sm = profile_sm / total
                ylabel = 'dP/dr  (cm⁻¹, normalised)'

            ax.plot(cy, profile_sm, color=color, linewidth=1.1,
                    label=f'{label}  (E={cx[ie]:.2g} eV)')

        ax.axvline(2.5, color='gray', linestyle=':', linewidth=0.8, label='Collimator 2.5 cm')
        ax.set_xlabel('Radial distance from beam axis (cm)')
        ax.set_ylabel(ylabel)
        ax.set_xlim(0, 3.0)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    _save(fig, 'ear2_spatial_radial_slices.png')


def plot_beam_size_vs_energy(vals_rb, cx, cy):
    mean_r  = np.zeros(len(cx))
    rms_r   = np.zeros(len(cx))
    r50     = np.zeros(len(cx))  # radius containing 50% of beam
    r90     = np.zeros(len(cx))  # radius containing 90% of beam

    for ie in range(len(cx)):
        profile = vals_rb[ie, :]
        total = profile.sum()
        if total <= 0:
            continue
        mean_r[ie] = np.sum(cy * profile) / total
        rms_r[ie]  = np.sqrt(np.sum(cy**2 * profile) / total)
        cumul = np.cumsum(profile) / total
        r50[ie] = cy[np.searchsorted(cumul, 0.50)]
        r90[ie] = cy[np.searchsorted(cumul, 0.90)]

    # Smooth curves over energy (they're noisy due to MC statistics)
    smooth = lambda x: uniform_filter1d(x, size=5)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    fig.suptitle('EAR2 Beam Spot Size vs Neutron Energy (FLUKA transport)', fontsize=13)

    ax1.plot(cx, smooth(mean_r), label='Mean r', color='tab:blue', linewidth=1.2)
    ax1.plot(cx, smooth(rms_r),  label='RMS r', color='tab:orange', linewidth=1.2)
    ax1.fill_between(cx, smooth(mean_r)-smooth(rms_r), smooth(mean_r)+smooth(rms_r),
                     alpha=0.15, color='tab:blue')
    ax1.axhline(2.5, color='gray', linestyle=':', linewidth=0.8, label='Collimator 2.5 cm')
    ax1.axvspan(E_X17_LOW_EV, E_X17_HIGH_EV, alpha=0.15, color='gold',
                label='X17 ROI 10–20 MeV')
    ax1.set_ylabel('Radius (cm)')
    ax1.set_ylim(0, 3.0)
    ax1.legend(fontsize=9)
    ax1.grid(True, which='both', alpha=0.3)
    _add_energy_ticks(ax1)

    ax2.plot(cx, smooth(r50), label='r₅₀  (50% enclosed)', color='tab:green', linewidth=1.2)
    ax2.plot(cx, smooth(r90), label='r₉₀  (90% enclosed)', color='tab:red', linewidth=1.2)
    ax2.axhline(2.5, color='gray', linestyle=':', linewidth=0.8, label='Collimator 2.5 cm')
    ax2.axvspan(E_X17_LOW_EV, E_X17_HIGH_EV, alpha=0.15, color='gold')
    ax2.set_xlabel('Neutron energy (eV)')
    ax2.set_ylabel('Radius (cm)')
    ax2.set_xscale('log')
    ax2.set_ylim(0, 3.0)
    ax2.legend(fontsize=9)
    ax2.grid(True, which='both', alpha=0.3)

    fig.tight_layout()
    _save(fig, 'ear2_spatial_beam_size.png')


def plot_target_fraction_vs_energy(vals_rb, cx, cy, ey):
    """Fraction of neutrons within radius r_target vs energy (uniform beam not assumed)."""
    fractions = {}
    for r_t in R_TARGETS_CM:
        ir = np.searchsorted(ey[1:], r_t)  # first bin edge > r_t
        frac = np.zeros(len(cx))
        for ie in range(len(cx)):
            profile = vals_rb[ie, :]
            total = profile.sum()
            if total > 0:
                frac[ie] = profile[:ir].sum() / total
        fractions[r_t] = frac

    smooth = lambda x: uniform_filter1d(x, size=5)

    fig, ax = plt.subplots(figsize=(11, 6))
    fig.suptitle('EAR2 Fraction of Beam Intercepted by Target vs Energy\n'
                 '(from FLUKA transport spatial profile)', fontsize=12)

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(R_TARGETS_CM)))
    for r_t, color in zip(R_TARGETS_CM, colors):
        frac = smooth(fractions[r_t])
        # Print X17 ROI average
        x17_mask = (cx >= E_X17_LOW_EV) & (cx <= E_X17_HIGH_EV)
        f_x17 = frac[x17_mask].mean() if x17_mask.any() else float('nan')
        ax.plot(cx, frac, color=color, linewidth=1.2,
                label=f'r < {r_t:.1f} cm   (X17 ROI avg: {f_x17*100:.0f}%)')

    ax.axvspan(E_X17_LOW_EV, E_X17_HIGH_EV, alpha=0.15, color='gold',
               label='X17 ROI 10–20 MeV')
    ax.set_xscale('log')
    ax.set_xlabel('Neutron energy (eV)')
    ax.set_ylabel('Fraction of beam within r_target')
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)
    _add_energy_ticks(ax)
    fig.tight_layout()
    _save(fig, 'ear2_spatial_target_fraction.png')


def _add_energy_ticks(ax):
    for val, label in [(1e-3, 'meV'), (1, 'eV'), (1e3, 'keV'), (1e6, 'MeV'), (1e8, '100 MeV')]:
        xlim = ax.get_xlim()
        if xlim[0] < val < xlim[1]:
            ax.axvline(val, color='white', linestyle=':', linewidth=0.4, alpha=0.4)


def _save(fig, name):
    path = OUT_DIR / name
    fig.savefig(path, dpi=150)
    print(f'Saved {path}')


if __name__ == '__main__':
    main()
