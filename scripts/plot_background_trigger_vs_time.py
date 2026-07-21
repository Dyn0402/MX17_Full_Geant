#!/usr/bin/env python3
"""
plot_background_trigger_vs_time.py — thermal-gate background vs arrival time
===========================================================================
Companion to the IPC in-gate spectrum.  Two panels, same arrival-time axis:

  (top)    the GENERAL background: neutron-capture rate vs time, by source
           (Al-dome nCapture, scintillator H-capture, total) — high statistics
           straight from the capture histograms.
  (bottom) what actually fires the singles (SiPM ∧ plastic) 2-arm coincidence,
           for several plastic thresholds b, with the IPC signal overlaid.

Why the shape is imported from the captures (as for the IPC reweighting): the
direct trigger MC has only ~40 in-gate background tags at b=0.5 (→0 at b=1.5),
but the capture γ energy is set by the nucleus, not the neutron arrival energy,
so the per-capture trigger efficiency is energy-independent.  We therefore take
the time SHAPE from the high-stats capture spectrum and the per-threshold
NORMALISATION from trigger_scan.json.  Physics refinement: H-capture (2.22 MeV
γ) deposits <0.5 MIP in the plastic, so for b ≥ 0.5 MIP the surviving
background is Al-dome-dominated (slightly later thermal peak).

Inputs (all local):
  docs/report/thermal_captures_subkev_full.npz   h_wall_al, h_scint_h(, edges)
  analysis/trigger_thermal/trigger_scan.json     leg2[a][b] background rates
  analysis/reweight/ipc_ingate_spectrum.npz      IPC signal reference
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NPZ  = ROOT / "docs/report/thermal_captures_subkev_full.npz"
TRIG = ROOT / "analysis/trigger_thermal/trigger_scan.json"
IPCF = ROOT / "analysis/reweight/ipc_ingate_spectrum.npz"
OUT  = ROOT / "analysis/reweight"

FLIGHT_MS = 1.41
GATE_MS   = 1.0
E_GATE    = (FLIGHT_MS / GATE_MS) ** 2
SIPM_A    = 0.5                              # SiPM leg threshold [MIP]
B_LINES   = [0.1, 0.3, 0.5, 0.7, 1.0]        # plastic thresholds [MIP]
H_KILL_B  = 0.5                              # b above which H-capture can't pass plastic
NFINE     = 1000

t_of_E = lambda E: FLIGHT_MS / np.sqrt(E)


def arrival_density(h, edges, w, Ef):
    """Convert a per-energy capture histogram to an arrival-time density on the
    shared fine energy grid Ef.  Returns (t_ms, dN/dt [per pulse per ms],
    in-gate total [per pulse]).  Straddling gate bin clipped to E<E_GATE."""
    E_lo = 10 ** edges[:-1]; E_hi = 10 ** edges[1:]; E_c = np.sqrt(E_lo * E_hi)
    frac = np.zeros_like(E_c)
    for i in range(len(E_c)):
        if E_lo[i] >= E_GATE:
            continue
        ehi = min(E_hi[i], E_GATE)
        frac[i] = np.log(ehi / E_lo[i]) / np.log(E_hi[i] / E_lo[i])
    total = float((h * frac).sum()) * w
    pop = (E_lo < E_GATE) & (h > 0)
    lnEc = np.log(E_c[pop]); dens = h[pop] / np.log(E_hi[pop] / E_lo[pop])
    densf = np.exp(np.interp(np.log(Ef), lnEc, np.log(dens),
                             left=-np.inf, right=-np.inf))   # 0 outside support
    tf = t_of_E(Ef)
    dNdt = densf * w * (2.0 / tf)                            # |dlnE/dt| = 2/t
    o = np.argsort(tf); tf, dNdt = tf[o], dNdt[o]
    integ = np.trapz(dNdt, tf)
    if integ > 0:
        dNdt *= total / integ
    return tf, dNdt, total


def main():
    d = np.load(NPZ)
    edges = d["loge_edges"]
    h_al  = d["h_wall_al"].astype(float)
    h_h   = d["h_scint_h"].astype(float)
    n_ev  = float(d["n_events"]); n_pp = float(d["n_per_pulse"]); w = n_pp / n_ev

    # shared fine energy grid over the populated in-gate range
    E_lo = 10 ** edges[:-1]
    pop  = (E_lo < E_GATE) & ((h_al + h_h) > 0)
    Emin = float(E_lo[pop].min())
    Ef   = np.exp(np.linspace(np.log(Emin), np.log(E_GATE), NFINE))

    t_al, al,  tot_al  = arrival_density(h_al,       edges, w, Ef)
    t_h,  hh,  tot_h   = arrival_density(h_h,        edges, w, Ef)
    t_t,  tot, tot_all = arrival_density(h_al + h_h, edges, w, Ef)
    # unit-area shapes for the trigger scaling
    shape_tot = tot / np.trapz(tot, t_t)
    shape_al  = al  / np.trapz(al,  t_al)

    # trigger normalisation per plastic threshold (2 full SiPM∧plastic legs)
    j = json.load(open(TRIG))
    a_grid = j["a_grid"]; b_grid = j["b_grid"]; Npt = j["n_pulse_thermal"]
    ai = a_grid.index(SIPM_A)
    leg2 = np.array(j["thermal"]["leg2"])
    def bg_total(b):
        bi = b_grid.index(round(b, 1))
        return leg2[ai][bi] * Npt, int(leg2[ai][bi] * n_ev)   # tags/pulse, MC evts

    # IPC signal reference
    z = np.load(IPCF)
    t_ipc = z["t_ms"]; ipc = z["dNdt_ipc_per_pulse_per_ms"]
    ipc_pp = float(z["ipc_per_pulse_ingate"])

    # plot range: trim the dead low-energy tail (density < 1e-3 of peak)
    t_max = float(t_t[tot > 1e-3 * tot.max()].max())
    msk = lambda t, y: (t <= t_max) & (y > 0)   # mask zeros for log axis

    # ── figure ────────────────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.6, 9.0), sharex=True)

    # (top) general background capture spectrum by source
    m = t_t <= t_max
    ax1.plot(t_t[m], tot[m], color="k", lw=2.2, label=f"total captures ({tot_all:.2e}/pulse)")
    ax1.plot(t_al[m], al[m], color="#1f77b4", lw=1.8, label=f"Al-dome nCapture ({tot_al:.2e}/pulse)")
    ax1.plot(t_h[m], hh[m], color="#d62728", lw=1.8, label=f"scint H-capture ({tot_h:.2e}/pulse)")
    ax1.fill_between(t_t[m], tot[m], color="0.85", zorder=0)
    ax1.set_ylabel("captures / pulse / ms")
    ax1.set_title("Thermal-gate background — general capture spectrum vs arrival time")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(fontsize=8.5, loc="upper right")
    ax1.set_ylim(bottom=0)

    # (bottom) what fires the singles SiPM∧plastic 2-arm trigger, per plastic thr.
    cmap = plt.cm.viridis(np.linspace(0.15, 0.85, len(B_LINES)))
    for b, col in zip(B_LINES, cmap):
        tp, ev = bg_total(b)
        shp = shape_al if b >= H_KILL_B else shape_tot   # H killed above ~0.5 MIP
        y = shp * tp; k = msk(t_t, y)
        ax2.plot(t_t[k], y[k], color=col, lw=1.9,
                 label=f"b={b:.1f} MIP:  {tp:.3f}/pulse  ({ev} MC evt)")
    k = msk(t_ipc, ipc)
    ax2.plot(t_ipc[k], ipc[k], color="limegreen", lw=2.4, ls="--",
             label=f"IPC signal ({ipc_pp:.2e}/pulse)")
    ax2.set_yscale("log")
    ax2.set_ylim(1e-8, 1.0)
    ax2.set_xlim(GATE_MS, t_max)
    ax2.set_xlabel("neutron arrival time  t [ms]   (gate: t > 1 ms)")
    ax2.set_ylabel("tags / pulse / ms")
    ax2.set_title(f"Background firing the singles (SiPM∧plastic) 2-arm trigger "
                  f"[SiPM ≥ {SIPM_A} MIP]  +  IPC signal")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(fontsize=8.0, loc="upper right", title="plastic threshold b "
               "(shape: Al-only for b≥0.5, Al+H below)")
    ax2.annotate("epithermal spill-in ≈ 0;  b ≥ 1.5 MIP → 0 tags (MC floor)",
                 (0.015, 0.03), xycoords="axes fraction", fontsize=8, color="0.4")

    fig.tight_layout()
    fig.savefig(OUT / "background_trigger_vs_time.pdf", bbox_inches="tight")
    fig.savefig(OUT / "background_trigger_vs_time.png", dpi=140, bbox_inches="tight")

    # ── small data dump (per-b totals + fine shapes) ──────────────────────────
    out = {
        "description": "Thermal-gate background vs arrival time: general capture "
                       "spectrum + singles(SiPM∧plastic) 2-arm trigger rate per "
                       "plastic threshold. Shape from high-stats captures, "
                       "normalisation from trigger_scan.json leg2.",
        "sipm_threshold_MIP": SIPM_A, "gate_ms": GATE_MS, "flight_ms": FLIGHT_MS,
        "epithermal_spillin": "≈0 (0 pair-tags in 5e8 epi neutrons)",
        "capture_per_pulse": {"al_dome": tot_al, "scint_H": tot_h, "total": tot_all},
        "ipc_signal_per_pulse": ipc_pp,
        "trigger_bg_per_pulse_vs_b": {
            f"{b:.1f}": {"tags_per_pulse": bg_total(b)[0], "mc_events": bg_total(b)[1],
                         "S_over_B": ipc_pp / bg_total(b)[0] if bg_total(b)[0] > 0 else None,
                         "shape": "Al-only" if b >= H_KILL_B else "Al+H"}
            for b in B_LINES},
        "note": "H-capture (2.22 MeV γ) deposits <0.5 MIP in plastic, so for "
                "b>=0.5 MIP the surviving background is Al-dome-dominated. The "
                "exact b-dependent shape/normalisation at high b is MC-floor "
                "limited (few tags) — confirm with the biased/high-stats run.",
    }
    (OUT / "background_trigger_vs_time.json").write_text(json.dumps(out, indent=1))

    print(f"captures/pulse in-gate: Al={tot_al:.3e}  H={tot_h:.3e}  total={tot_all:.3e}")
    print(f"IPC signal/pulse      : {ipc_pp:.3e}")
    print("\nbackground trigger (SiPM≥%.1f) vs plastic threshold b:" % SIPM_A)
    for b in B_LINES:
        tp, ev = bg_total(b)
        sb = ipc_pp / tp if tp > 0 else float('inf')
        print(f"  b={b:.1f} MIP: {tp:.4f} tags/pulse ({ev:4d} MC evt)   S/B = {sb:.2e}")
    print(f"\nwrote -> {OUT}/background_trigger_vs_time.pdf/.png/.json")


if __name__ == "__main__":
    main()
