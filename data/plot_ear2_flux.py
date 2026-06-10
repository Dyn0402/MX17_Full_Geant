#!/usr/bin/env python3
"""
Plot EAR2 neutron flux profiles from the evaluated n_TOF Phase 3 flux files.
Data from: /eos/experiment/ntof/flux_RF_BIF/FLUX/EAR2/

Plots:
  1. All four flux unit variants (log-log)
  2. Ph3 evaluated flux vs older 2014-2015 evaluation
  3. TOF-space histogram
  4. Collimator radial footprint context (r=2.5 cm assumed for cm² conversion)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import uproot
from pathlib import Path

DATA_DIR = Path('/home/dylan/CLionProjects/MX17_Full_Geant/data')
OUT_DIR = Path('/home/dylan/CLionProjects/MX17_Full_Geant/data')

# X17 search energy range of interest
E_X17_LOW_MEV = 10.0   # MeV
E_X17_HIGH_MEV = 20.0  # MeV

# EAR2 collimator radius used for neutrons/cm² conversion
R_COLLIMATOR_CM = 2.5
AREA_CM2 = np.pi * R_COLLIMATOR_CM**2


def load_hist(f, key):
    h = f[key]
    vals = h.values()
    edges = h.axis(0).edges()
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, edges, vals


def main():
    plot_flux_units()
    plot_ph3_vs_old()
    plot_tof()
    plot_radial_context()
    plt.show()


def plot_flux_units():
    with uproot.open(DATA_DIR / 'fluxEAR2-Ph3_in_different_units.root') as f:
        E_isol, _, flux_isol = load_hist(f, 'flux_n_pulse_isolet_100bpd;1')
        E_noiso, _, flux_noiso = load_hist(f, 'flux_n_pulse_NOisolet_100bpd;1')
        E_cm2p, _, flux_cm2p = load_hist(f, 'flux_n_cm2_pulse_100bpd;1')
        E_cm2s, _, flux_cm2s = load_hist(f, 'flux_n_cm2_s_100bpd;1')

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle('EAR2 Neutron Flux — n_TOF Phase 3 Evaluated', fontsize=14)

    panels = [
        (axes[0, 0], E_isol,  flux_isol,  'n / eV / pulse',    'tab:blue'),
        (axes[0, 1], E_noiso, flux_noiso,  'n / pulse',          'tab:orange'),
        (axes[1, 0], E_cm2p,  flux_cm2p,  'n / cm² / pulse',   'tab:green'),
        (axes[1, 1], E_cm2s,  flux_cm2s,  'n / cm² / s',        'tab:red'),
    ]

    for ax, E, flux, ylabel, color in panels:
        ax.plot(E, flux, color=color, linewidth=0.8)
        ax.axvspan(E_X17_LOW_MEV * 1e6, E_X17_HIGH_MEV * 1e6,
                   alpha=0.15, color='gold', label=f'X17 ROI {E_X17_LOW_MEV}–{E_X17_HIGH_MEV} MeV')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Neutron energy (eV)')
        ax.set_ylabel(ylabel)
        ax.set_xlim(1e-3, 1e8)
        ax.legend(fontsize=8)
        ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
        _add_energy_labels(ax)
        ax.grid(True, which='both', alpha=0.3)

    fig.tight_layout()
    out = OUT_DIR / 'ear2_flux_units.png'
    fig.savefig(out, dpi=150)
    print(f'Saved {out}')


def plot_ph3_vs_old():
    with uproot.open(DATA_DIR / 'fluxEAR2-Ph3_in_different_units.root') as f:
        E_ph3, _, flux_ph3 = load_hist(f, 'flux_n_pulse_isolet_100bpd;1')

    with uproot.open(DATA_DIR / 'nTOF-Ph3_fluence_EAR2_2014-2015.root') as f:
        E_old, _, flux_old = load_hist(f, 'hFlux_eval_ear2_2014_2015_100bpd;1')

    fig, (ax_main, ax_ratio) = plt.subplots(2, 1, figsize=(11, 8),
                                             gridspec_kw={'height_ratios': [3, 1.5]},
                                             sharex=True)
    fig.suptitle('EAR2 Isolethargic Flux: Ph3 evaluated vs 2014–2015', fontsize=13)

    ax_main.plot(E_ph3, flux_ph3, color='tab:blue', linewidth=0.9, label='Ph3 evaluated (current)')
    ax_main.plot(E_old, flux_old, color='tab:orange', linewidth=0.9,
                 linestyle='--', alpha=0.8, label='2014–2015 commissioning')
    ax_main.axvspan(E_X17_LOW_MEV * 1e6, E_X17_HIGH_MEV * 1e6,
                    alpha=0.15, color='gold', label=f'X17 ROI {E_X17_LOW_MEV}–{E_X17_HIGH_MEV} MeV')
    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_ylabel('n / eV / pulse')
    ax_main.legend()
    ax_main.grid(True, which='both', alpha=0.3)
    _add_energy_labels(ax_main)

    # Ratio: interpolate old onto Ph3 grid
    ratio = np.where(flux_old > 0, flux_ph3 / flux_old, np.nan)
    ax_ratio.plot(E_ph3, ratio, color='black', linewidth=0.8)
    ax_ratio.axhline(1.0, color='gray', linestyle='--', linewidth=0.8)
    ax_ratio.axvspan(E_X17_LOW_MEV * 1e6, E_X17_HIGH_MEV * 1e6, alpha=0.15, color='gold')
    ax_ratio.set_xscale('log')
    ax_ratio.set_xlabel('Neutron energy (eV)')
    ax_ratio.set_ylabel('Ph3 / 2014–2015')
    ax_ratio.set_ylim(0.5, 2.5)
    ax_ratio.grid(True, which='both', alpha=0.3)

    fig.tight_layout()
    out = OUT_DIR / 'ear2_flux_ph3_vs_old.png'
    fig.savefig(out, dpi=150)
    print(f'Saved {out}')


def plot_tof():
    with uproot.open(DATA_DIR / 'nTOF-Ph3_fluence_EAR2_2014-2015_TOF-E.root') as f:
        E_hires, _, flux_hires = load_hist(f, 'hE_eval_ear2_2014_2015_1000bpd;1')
        tof_centers, tof_edges, flux_tof = load_hist(f, 'hTOF_SiMon2pos_ear2_2014_2015_1000bpd;1')

    # TOF axis is in eV (edges are in eV, representing TOF-space binned energy)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle('EAR2 Flux — TOF-E file (1000 bins/decade, 2014–2015)', fontsize=13)

    ax = axes[0]
    ax.plot(E_hires, flux_hires, color='tab:purple', linewidth=0.6)
    ax.axvspan(E_X17_LOW_MEV * 1e6, E_X17_HIGH_MEV * 1e6,
               alpha=0.15, color='gold', label=f'X17 ROI')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Neutron energy (eV)')
    ax.set_ylabel('n / eV / pulse')
    ax.set_title('High-resolution energy flux (1000 bpd)')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    _add_energy_labels(ax)

    ax = axes[1]
    ax.plot(tof_centers, flux_tof, color='tab:brown', linewidth=0.6)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Neutron energy (eV)  [TOF binning]')
    ax.set_ylabel('n / eV / pulse')
    ax.set_title('TOF-space binned flux (SiMon2 position)')
    ax.grid(True, which='both', alpha=0.3)
    _add_energy_labels(ax)

    fig.tight_layout()
    out = OUT_DIR / 'ear2_flux_tof.png'
    fig.savefig(out, dpi=150)
    print(f'Saved {out}')


def plot_radial_context():
    """
    Show integrated flux through collimator vs energy, and the implied
    'effective beam radius' (fixed at 2.5 cm per the Ph3 evaluation assumptions).
    Also shows the fraction of the full flux intercepted by a smaller target of
    varying radius, assuming a uniform beam profile within the collimator.
    """
    with uproot.open(DATA_DIR / 'fluxEAR2-Ph3_in_different_units.root') as f:
        E, _, flux_per_pulse = load_hist(f, 'flux_n_pulse_NOisolet_100bpd;1')
        _, _, flux_cm2 = load_hist(f, 'flux_n_cm2_pulse_100bpd;1')

    # Integrated flux above each energy threshold
    cumulative = np.flip(np.cumsum(np.flip(flux_per_pulse)))

    # Flux density (n/cm²/pulse) reconstructed from ratio — should be flat = 1/area
    # flux_cm2 = flux_per_pulse / area  (uniform beam within r_coll assumed)
    # We can show neutrons intercepted by a target of radius r < r_coll (uniform beam)
    r_vals_cm = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
    fractions = (r_vals_cm / R_COLLIMATOR_CM)**2  # uniform beam: fraction = (r/r_coll)²

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle('EAR2 Radial / Collimator Context (r_coll = 2.5 cm)', fontsize=13)

    # Left: cumulative flux above E
    ax = axes[0]
    ax.plot(E, cumulative, color='tab:blue', linewidth=0.9)
    ax.axvspan(E_X17_LOW_MEV * 1e6, E_X17_HIGH_MEV * 1e6,
               alpha=0.2, color='gold', label=f'X17 ROI {E_X17_LOW_MEV}–{E_X17_HIGH_MEV} MeV')
    x17_mask = (E >= E_X17_LOW_MEV * 1e6) & (E <= E_X17_HIGH_MEV * 1e6)
    n_x17 = float(np.sum(flux_per_pulse[x17_mask]))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Neutron energy (eV)')
    ax.set_ylabel('n / pulse  [E > threshold]')
    ax.set_title('Cumulative flux above energy threshold\n(full 2.5 cm collimator)')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    ax.text(0.05, 0.05,
            f'X17 ROI total: {n_x17:.2e} n/pulse\n(r < {R_COLLIMATOR_CM} cm)',
            transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.6))
    _add_energy_labels(ax)

    # Right: neutrons intercepted by smaller target (uniform beam)
    ax = axes[1]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(r_vals_cm)))
    for r, frac, color in zip(r_vals_cm, fractions, colors):
        ax.plot(E, flux_per_pulse * frac, color=color, linewidth=0.9,
                label=f'r = {r:.1f} cm  ({frac*100:.0f}%)')
    ax.axvspan(E_X17_LOW_MEV * 1e6, E_X17_HIGH_MEV * 1e6, alpha=0.15, color='gold')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Neutron energy (eV)')
    ax.set_ylabel('n / pulse / bin')
    ax.set_title('Flux intercepted by target of radius r\n(uniform beam within collimator assumed)')
    ax.legend(fontsize=8)
    ax.grid(True, which='both', alpha=0.3)
    _add_energy_labels(ax)

    fig.tight_layout()
    out = OUT_DIR / 'ear2_flux_radial_context.png'
    fig.savefig(out, dpi=150)
    print(f'Saved {out}')
    print(f'\nNote: The "radial footprint vs energy" (beam profile variation with E) is not in these files.')
    print(f'That information lives in the ntof/transport simulation outputs.')
    print(f'These plots assume a uniform beam within r=2.5 cm at all energies.')


def _add_energy_labels(ax):
    xlim = ax.get_xlim()
    labels = {'meV': 1e-3, 'eV': 1, 'keV': 1e3, 'MeV': 1e6}
    ymin, ymax = ax.get_ylim()
    for label, val in labels.items():
        if xlim[0] < val < xlim[1]:
            ax.axvline(val, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)


if __name__ == '__main__':
    main()
