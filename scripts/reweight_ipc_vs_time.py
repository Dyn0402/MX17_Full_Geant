#!/usr/bin/env python3
"""
reweight_ipc_vs_time.py — high-statistics IPC production vs arrival time
========================================================================
Method #2 in the thermal handoff discussion (2026-07-21): the radiative
capture rate is a fixed analytic multiple of the abundant (n,p)t absorption,

    R_nγ(E) = R_np(E) × σ_nγ(E)/σ_np(E)         [ENDF/B-VIII.0]

so the *shape* of IPC production vs neutron energy (= arrival time via TOF)
is known from the (n,p)t sample — 7.4×10⁸ events — with negligible
statistical error, instead of the 15 raw direct captures.

Inputs (all local, no EOS):
  docs/report/thermal_captures_subkev_full.npz   h_gas_np(E), h_gas_rad(E)
  data/He3.h5                                    σ_nγ (MT=102), σ_np (MT=103)

Arrival time (EAR2 → target, 19.5 m):  t[ms] = 1.41 / sqrt(E[eV])
  gate t > 1 ms  ⇔  E_n < 1.99 eV.

Outputs:
  analysis/reweight/ipc_vs_time.json   numbers for operational decisions
  analysis/reweight/ipc_vs_time.pdf    differential + cumulative figure
"""
import json
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NPZ  = ROOT / "docs/report/thermal_captures_subkev_full.npz"
H5   = ROOT / "data/He3.h5"
OUT  = ROOT / "analysis/reweight"
OUT.mkdir(parents=True, exist_ok=True)

# ── beam / detector constants (match analyze_thermal_captures.py) ──────────────
FLIGHT_MS   = 1.41          # t[ms] = FLIGHT_MS / sqrt(E[eV]) at 19.5 m
GATE_MS     = 1.0
PULSES_DAY  = 1.929e4       # EAR2 pulses/day (7e12 ppp)
SIGMA_RATIO_TH = 54e-6 / 5333.0     # flat thermal σ_nγ/σ_np = 1.013e-8 (floor)


def load_xs(mt):
    """Return (E[eV], xs[barn]) for ENDF reaction MT in the OpenMC He3.h5.
    Structure (verified): energy grid He3/energy/294K [eV]; per-reaction xs at
    He3/reactions/reaction_<MT>/294K/xs, sliced by attr 'threshold_idx'."""
    with h5py.File(H5, "r") as f:
        g = f["He3"]
        E_eV = g["energy/294K"][...]                 # already in eV
        ds   = g[f"reactions/reaction_{mt}/294K/xs"]
        xs   = ds[...]
        i0   = ds.attrs.get("threshold_idx", 0)
    return np.asarray(E_eV[i0:i0 + len(xs)]), np.asarray(xs)


def loglog_interp(xq, x, y):
    lx, ly = np.log(x), np.log(np.clip(y, 1e-300, None))
    return np.exp(np.interp(np.log(xq), lx, ly))


def main():
    d = np.load(NPZ)
    edges = d["loge_edges"]                       # log10(E/eV), 61 edges
    h_np  = d["h_gas_np"].astype(float)           # (n,p)t counts / bin  (abundant)
    h_rad = d["h_gas_rad"].astype(float)          # direct (n,γ) counts / bin (15)
    n_ev  = float(d["n_events"])                  # neutrons simulated (1e9)
    n_pp  = float(d["n_per_pulse"])               # beam n/pulse in this window
    a_ipc = float(d["alpha_ipc"])                 # IPC pairs / radiative capture
    br    = float(d["br_x17"])                    # X17 / IPC pair
    w     = n_pp / n_ev                           # sim-count → per-pulse

    Ec = 10 ** (0.5 * (edges[:-1] + edges[1:]))   # bin-centre energy [eV]
    t  = FLIGHT_MS / np.sqrt(Ec)                   # arrival time [ms]

    # σ_nγ/σ_np(E) evaluated at each bin centre (loglog ENDF interpolation)
    E_np, xs_np = load_xs("103")
    E_ng, xs_ng = load_xs("102")
    ratio = loglog_interp(Ec, E_ng, xs_ng) / loglog_interp(Ec, E_np, xs_np)

    # ── the reweighting ───────────────────────────────────────────────────────
    rad_rw   = h_np * ratio                         # expected radiative capt / bin
    cap_pp   = rad_rw * w                           # captures / pulse / bin
    ipc_pp   = cap_pp * a_ipc                        # IPC pairs / pulse / bin
    x17_pp   = ipc_pp * br                           # X17 / pulse / bin

    # flat-floor comparison (single thermal ratio, the old "floor" method)
    ipc_pp_flat = h_np * SIGMA_RATIO_TH * w * a_ipc

    in_gate = t > GATE_MS                             # E < 1.99 eV
    tot_ipc      = ipc_pp.sum()
    gate_ipc     = ipc_pp[in_gate].sum()
    gate_ipc_flat = ipc_pp_flat[in_gate].sum()
    gate_cap     = cap_pp[in_gate].sum()

    # effective statistics: Poisson error of the reweighted estimate is driven
    # by the (n,p)t counts, not the 15 direct captures.
    eff_N_gate = h_np[in_gate].sum()
    frac_err_gate = 1.0 / np.sqrt(eff_N_gate) if eff_N_gate else np.inf

    # cumulative IPC/pulse arriving LATER than t (i.e. gate opens at t) —
    # the operational curve: integrate from high-E(early) to low-E(late).
    order = np.argsort(t)
    t_sorted   = t[order]
    ipc_sorted = ipc_pp[order]
    cum_late   = np.cumsum(ipc_sorted[::-1])[::-1]   # sum of all bins with t' >= t

    # ── report ────────────────────────────────────────────────────────────────
    print(f"(n,p)t events (all)           : {h_np.sum():.3e}")
    print(f"(n,p)t events in-gate (>1 ms) : {eff_N_gate:.3e}  -> stat err {frac_err_gate*100:.2f}%")
    print(f"direct (n,γ) observed         : {int(h_rad.sum())} total, {int(h_rad[in_gate].sum())} in-gate")
    print(f"per-pulse scale w             : {w:.3e}")
    print("-" * 64)
    print(f"radiative captures / pulse    : total {cap_pp.sum():.3e}  in-gate {gate_cap:.3e}")
    print(f"IPC pairs / pulse (E-dep)     : total {tot_ipc:.3e}  in-gate {gate_ipc:.3e}")
    print(f"IPC pairs / pulse (flat floor): total {ipc_pp_flat.sum():.3e}  in-gate {gate_ipc_flat:.3e}")
    print(f"  -> E-dep / flat (in-gate)   : {gate_ipc/gate_ipc_flat:.2f}x")
    print(f"IPC / day  (in-gate)          : {gate_ipc*PULSES_DAY:.2f}  +/- {gate_ipc*PULSES_DAY*frac_err_gate:.2f} (stat)")
    print(f"X17 produced / day (in-gate)  : {gate_ipc*br*PULSES_DAY:.3f}")
    # cross-check: reweighted total captures vs the 15 observed directly
    pred_capt_sim = rad_rw.sum()
    print("-" * 64)
    print(f"CROSS-CHECK reweighted predicted direct captures in 1e9 sim: "
          f"{pred_capt_sim:.1f}  vs {int(h_rad.sum())} observed "
          f"(Poisson {h_rad.sum():.0f}±{np.sqrt(h_rad.sum()):.1f})")

    out = {
        "method": "reweight (n,p)t x sigma_ng/sigma_np(E)",
        "gate_ms": GATE_MS, "flight_ms": FLIGHT_MS, "pulses_per_day": PULSES_DAY,
        "eff_N_np_in_gate": float(eff_N_gate),
        "stat_frac_err_in_gate": float(frac_err_gate),
        "direct_obs_total": int(h_rad.sum()), "direct_obs_in_gate": int(h_rad[in_gate].sum()),
        "cap_per_pulse_in_gate": float(gate_cap),
        "ipc_per_pulse_in_gate": float(gate_ipc),
        "ipc_per_pulse_in_gate_flatfloor": float(gate_ipc_flat),
        "ipc_per_day_in_gate": float(gate_ipc * PULSES_DAY),
        "x17_per_day_in_gate": float(gate_ipc * br * PULSES_DAY),
        "reweight_predicted_direct_capt_in_sim": float(pred_capt_sim),
        "alpha_ipc": a_ipc, "br_x17": br,
        "t_ms": t.tolist(), "E_eV": Ec.tolist(),
        "ipc_per_pulse_bin": ipc_pp.tolist(),
        "ipc_per_pulse_bin_flat": ipc_pp_flat.tolist(),
        "cum_ipc_per_pulse_later_than_t": {"t_ms": t_sorted.tolist(),
                                           "cum_ipc": cum_late.tolist()},
    }
    (OUT / "ipc_vs_time.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {OUT/'ipc_vs_time.json'}")

    # ── figure ────────────────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 8.2))

    # top: differential IPC production per pulse vs arrival time
    ax1.step(t, ipc_pp, where="mid", color="#1f77b4", lw=1.8,
             label=r"reweighted, $\sigma_{n\gamma}/\sigma_{np}(E)$ (7.4$\times10^8$ eff.)")
    ax1.step(t, ipc_pp_flat, where="mid", color="#888", lw=1.2, ls="--",
             label=r"flat thermal ratio (old floor)")
    # direct observed events, scaled to per-pulse-per-bin, as sparse crosses
    rad_pp_obs = h_rad * w * a_ipc
    m = rad_pp_obs > 0
    ax1.plot(t[m], rad_pp_obs[m], "rx", ms=8, mew=2,
             label=f"direct obs ({int(h_rad.sum())} events)")
    ax1.axvspan(GATE_MS, t.max() * 1.5, color="0.9", zorder=0)
    ax1.axvline(GATE_MS, color="k", ls=":", lw=1)
    ax1.text(GATE_MS * 1.05, ax1.get_ylim()[1], " gate: t>1 ms\n (E<1.99 eV)",
             va="top", fontsize=8.5)
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("arrival time  t [ms]   (= 1.41/√E,  19.5 m)")
    ax1.set_ylabel("IPC pairs / pulse / bin")
    ax1.set_title("Thermal-gate IPC production vs arrival time (reweighted)")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(fontsize=8.5, loc="lower left")
    secx = ax1.secondary_xaxis(
        "top", functions=(lambda tt: (FLIGHT_MS / np.clip(tt, 1e-9, None))**2,
                          lambda E: FLIGHT_MS / np.sqrt(np.clip(E, 1e-30, None))))
    secx.set_xlabel("neutron energy  E [eV]")

    # bottom: cumulative IPC/pulse for a gate opening at t (operational curve)
    ax2.step(t_sorted, cum_late, where="post", color="#d62728", lw=2)
    ax2.axvline(GATE_MS, color="k", ls=":", lw=1)
    ax2.plot([GATE_MS], [gate_ipc], "ko", ms=7)
    ax2.annotate(f"  {gate_ipc:.2e} IPC/pulse\n  ({gate_ipc*PULSES_DAY:.1f} IPC/day)",
                 (GATE_MS, gate_ipc), fontsize=9, va="center")
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel("gate-open time  t [ms]")
    ax2.set_ylabel("IPC pairs / pulse arriving after t")
    ax2.set_title("Cumulative in-gate IPC yield vs gate-open time")
    ax2.grid(True, which="both", alpha=0.25)

    fig.tight_layout()
    fig.savefig(OUT / "ipc_vs_time.pdf", bbox_inches="tight")
    fig.savefig(OUT / "ipc_vs_time.png", dpi=140, bbox_inches="tight")
    print(f"wrote {OUT/'ipc_vs_time.pdf'}")


if __name__ == "__main__":
    main()
