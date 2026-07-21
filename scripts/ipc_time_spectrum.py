#!/usr/bin/env python3
"""
ipc_time_spectrum.py — loader/sampler for the in-gate IPC arrival-time spectrum
===============================================================================
Reusable API over the data package written by export_ipc_ingate_spectrum.py.
Depends only on numpy.  Use it to plot the spectrum, evaluate the production
density at any time, or generate arbitrarily many *unbinned* arrival times /
neutron energies distributed exactly as the reweighted spectrum.

    from ipc_time_spectrum import IPCTimeSpectrum
    s = IPCTimeSpectrum()                      # loads analysis/reweight/*.npz
    s.ipc_per_pulse, s.ipc_per_day             # absolute in-gate normalisation
    y  = s.density(t_ms)                        # IPC/pulse/ms at arbitrary t
    F  = s.cdf_at(t_ms)                         # fraction of in-gate IPC with arr<=t
    t  = s.sample(100000)                       # unbinned arrival times [ms]
    E  = s.sample_energy(100000)               # matching neutron energies [eV]
    n_expected = s.ipc_per_pulse * n_pulses     # scale to a run

Rescale to captures or X17:  captures = IPC / s.alpha_ipc ;  X17 = IPC * s.br_x17.
"""
from pathlib import Path
import json

import numpy as np

_DEFAULT = Path(__file__).resolve().parent.parent / "analysis/reweight/ipc_ingate_spectrum.npz"


class IPCTimeSpectrum:
    def __init__(self, npz_path=_DEFAULT):
        z = np.load(npz_path, allow_pickle=False)
        self.t_ms  = z["t_ms"]                              # ascending
        self.dNdt  = z["dNdt_ipc_per_pulse_per_ms"]         # IPC/pulse/ms
        self.cdf   = z["cdf"]                               # 0..1
        self.E_eV  = z["E_eV"]                              # matched to t_ms
        self.ipc_per_pulse = float(z["ipc_per_pulse_ingate"])
        self.ipc_per_day   = float(z["ipc_per_day_ingate"])
        self.captures_per_pulse = float(z["captures_per_pulse_ingate"])
        self.alpha_ipc = float(z["alpha_ipc"])
        self.br_x17    = float(z["br_x17"])
        self.pulses_per_day = float(z["pulses_per_day"])
        self.flight_ms = float(z["flight_ms"])
        self.gate_ms   = float(z["gate_ms"])
        self.meta = json.loads(str(z["meta_json"]))
        # raw bins (for exact folding)
        self.bin_t_lo = z["bin_t_lo"]; self.bin_t_hi = z["bin_t_hi"]
        self.bin_ipc_per_pulse = z["bin_ipc_per_pulse"]

    # ── continuous evaluation ────────────────────────────────────────────────
    def density(self, t_ms):
        """IPC pairs / pulse / ms at arrival time(s) t_ms (0 outside the gate)."""
        return np.interp(t_ms, self.t_ms, self.dNdt, left=0.0, right=0.0)

    def cdf_at(self, t_ms):
        """Fraction of in-gate IPC production with arrival time <= t_ms."""
        return np.interp(t_ms, self.t_ms, self.cdf, left=0.0, right=1.0)

    def E_at(self, t_ms):
        """Neutron energy [eV] for a given arrival time [ms]."""
        return (self.flight_ms / np.asarray(t_ms)) ** 2

    # ── unbinned sampling (inverse-transform on the fine CDF) ─────────────────
    def sample(self, n, rng=None):
        """n unbinned arrival times [ms] ~ the in-gate production spectrum."""
        rng = rng or np.random.default_rng()
        return np.interp(rng.random(int(n)), self.cdf, self.t_ms)

    def sample_energy(self, n, rng=None):
        """n unbinned neutron energies [eV] ~ the in-gate production spectrum."""
        return (self.flight_ms / self.sample(n, rng)) ** 2

    def expected_counts(self, n_pulses):
        """Expected number of produced in-gate IPC pairs for n_pulses."""
        return self.ipc_per_pulse * n_pulses


if __name__ == "__main__":
    s = IPCTimeSpectrum()
    print("In-gate IPC production spectrum")
    print(f"  ∫ density        = {s.ipc_per_pulse:.4e} IPC/pulse "
          f"({s.ipc_per_day:.2f} IPC/day)")
    print(f"  time span        = {s.t_ms.min():.2f} – {s.t_ms.max():.2f} ms")
    print(f"  peak at          = {s.t_ms[np.argmax(s.dNdt)]:.2f} ms")
    # sampling self-check: sampled mean time vs density-weighted mean
    t = s.sample(200000)
    tbar_num = np.trapz(s.t_ms * s.dNdt, s.t_ms) / np.trapz(s.dNdt, s.t_ms)
    print(f"  <t> sampled      = {t.mean():.3f} ms   (curve {tbar_num:.3f} ms)")
    print(f"  median t sampled = {np.median(t):.3f} ms   "
          f"(curve {np.interp(0.5, s.cdf, s.t_ms):.3f} ms)")
