#!/usr/bin/env python3
"""
export_ipc_ingate_spectrum.py — in-gate (t > 1 ms) IPC production spectrum
==========================================================================
Companion to reweight_ipc_vs_time.py.  Restricts the reweighted ³He(n,γ)→IPC
production to the thermal gate (arrival time t > 1 ms ⇔ E_n < 1.99 eV), draws
it in LINEAR-LINEAR scale, and — the main point — writes a self-describing,
near-unbinned data package other analyses can plot or fold in:

  ipc_ingate_spectrum.npz        everything (fine grid + bins + scalars + meta)
  ipc_ingate_spectrum.json       portable metadata + tables (no numpy needed)
  ipc_ingate_spectrum_fine.csv   t_ms, dN/dt [IPC/pulse/ms], CDF, E_eV
                                 (plot directly, or inverse-CDF sample)
  ipc_ingate_spectrum_bins.csv   raw reweighted histogram + Poisson errors
  ipc_vs_time_ingate_linlin.pdf/png

"Unbinned" strategy: the reweighting's intrinsic resolution is the 10-bins/
decade (n,p)t histogram, but the spectrum is smooth (1/v cross-sections, no
resonances), so we ship a fine (1200-pt) differential density + its normalised
CDF.  Downstream can evaluate the density at any t, or inverse-transform the
CDF to generate arbitrarily many unbinned arrival times/energies — see
scripts/ipc_time_spectrum.py for the loader/sampler.
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

# ── constants (match reweight_ipc_vs_time.py / analyze_thermal_captures.py) ────
FLIGHT_MS  = 1.41            # t[ms] = FLIGHT_MS / sqrt(E[eV]) at 19.5 m
GATE_MS    = 1.0
PULSES_DAY = 1.929e4
E_GATE     = (FLIGHT_MS / GATE_MS) ** 2      # 1.988 eV
NFINE      = 1200
DATE       = "2026-07-21"
PROVENANCE = ("reweight of docs/report/thermal_captures_subkev_full.npz "
              "(Run B sub-keV, 1e9 EAR2 neutrons, 2026-06-11) × ENDF/B-VIII.0 "
              "σ_nγ/σ_np(E) from data/He3.h5; validated against a nCapture-"
              "biased Geant4 run (7–8e-9 capt/n vs 7.3e-9 reweighted).")

t_of_E = lambda E: FLIGHT_MS / np.sqrt(E)
E_of_t = lambda t: (FLIGHT_MS / t) ** 2


def load_xs(mt):
    with h5py.File(H5, "r") as f:
        g  = f["He3"]
        E  = g["energy/294K"][...]
        ds = g[f"reactions/reaction_{mt}/294K/xs"]
        xs = ds[...]
        i0 = ds.attrs.get("threshold_idx", 0)
    return np.asarray(E[i0:i0 + len(xs)]), np.asarray(xs)


def loglog(xq, x, y):
    return np.exp(np.interp(np.log(xq), np.log(x), np.log(np.clip(y, 1e-300, None))))


_E103, _X103 = load_xs("103")   # (n,p)
_E102, _X102 = load_xs("102")   # (n,γ)
def ratio(E):
    return loglog(E, _E102, _X102) / loglog(E, _E103, _X103)


def main(npz_path=NPZ, out_dir=OUT, geom_tag=""):
    global OUT
    OUT = out_dir
    OUT.mkdir(parents=True, exist_ok=True)
    d     = np.load(npz_path)
    edges = d["loge_edges"]
    h_np  = d["h_gas_np"].astype(float)
    n_ev  = float(d["n_events"]); n_pp = float(d["n_per_pulse"])
    a_ipc = float(d["alpha_ipc"]); br = float(d["br_x17"])
    w     = n_pp / n_ev

    E_lo = 10 ** edges[:-1]; E_hi = 10 ** edges[1:]
    E_c  = np.sqrt(E_lo * E_hi)

    # ── per-bin reweighted IPC/pulse (clip the gate-straddling bin to E<E_GATE),
    #    using the log-uniform-in-E bin-averaged ENDF ratio for accuracy ────────
    NS = 64
    ipc_bin = np.zeros_like(E_c)
    rat_bin = np.zeros_like(E_c)
    frac_ig = np.zeros_like(E_c)
    for i in range(len(E_c)):
        if E_lo[i] >= E_GATE:
            continue
        ehi   = min(E_hi[i], E_GATE)
        frac  = np.log(ehi / E_lo[i]) / np.log(E_hi[i] / E_lo[i])   # in-gate frac
        Es    = np.exp(np.linspace(np.log(E_lo[i]), np.log(ehi), NS))
        rbar  = ratio(Es).mean()
        rat_bin[i] = rbar; frac_ig[i] = frac
        ipc_bin[i] = h_np[i] * frac * rbar * w * a_ipc

    ingate      = E_lo < E_GATE
    ipc_ingate  = float(ipc_bin[ingate].sum())          # ground-truth normalisation
    cap_ingate  = ipc_ingate / a_ipc
    x17_ingate  = ipc_ingate * br
    # effective statistics behind the shape (Poisson on the abundant (n,p)t)
    Nnp_ingate  = float((h_np * frac_ig)[ingate].sum())

    # ── fine differential density for smooth plot + sampling ───────────────────
    pop  = ingate & (h_np > 0)
    Emin = float(E_lo[pop].min())
    Ef   = np.exp(np.linspace(np.log(Emin), np.log(E_GATE), NFINE))
    lnEc = np.log(E_c[pop])
    dens_np = h_np[pop] / np.log(E_hi[pop] / E_lo[pop])      # (n,p)t per unit lnE
    dens_np_f = np.exp(np.interp(np.log(Ef), lnEc, np.log(dens_np)))
    dNdlnE  = dens_np_f * ratio(Ef) * w * a_ipc              # IPC/pulse per lnE
    tf      = t_of_E(Ef)
    dNdt    = dNdlnE * (2.0 / tf)                            # IPC/pulse per ms

    o  = np.argsort(tf)                                      # ascending time
    tf, dNdt, Ef = tf[o], dNdt[o], Ef[o]
    # renormalise the smooth curve so its integral equals the exact bin sum
    dNdt *= ipc_ingate / np.trapz(dNdt, tf)
    cdf = np.concatenate([[0.0],
          np.cumsum(0.5 * (dNdt[1:] + dNdt[:-1]) * np.diff(tf))])
    cdf /= cdf[-1]

    # ── scalars / metadata ────────────────────────────────────────────────────
    meta = {
        "description": "In-gate (t>1 ms) IPC production vs neutron arrival time",
        "date": DATE, "provenance": PROVENANCE,
        "time_energy_relation": "t_ms = FLIGHT_MS / sqrt(E_eV); E_eV = (FLIGHT_MS/t_ms)**2",
        "FLIGHT_MS": FLIGHT_MS, "flight_length_m": 19.5,
        "gate_ms": GATE_MS, "E_gate_eV": E_GATE,
        "pulses_per_day": PULSES_DAY,
        "alpha_ipc_per_capture": a_ipc, "br_x17_per_ipc": br,
        "units": {"dNdt": "IPC pairs / pulse / ms", "cdf": "fraction of in-gate IPC with arrival time <= t",
                  "ipc_per_pulse": "IPC pairs per beam pulse", "bin ipc": "IPC pairs / pulse (integrated over bin)"},
        "normalisation": {
            "ipc_per_pulse_ingate": ipc_ingate,
            "ipc_per_day_ingate":   ipc_ingate * PULSES_DAY,
            "captures_per_pulse_ingate": cap_ingate,
            "x17_per_day_ingate":   x17_ingate * PULSES_DAY,
            "eff_np_counts_ingate": Nnp_ingate,
            "shape_stat_frac_err":  1.0 / np.sqrt(Nnp_ingate),
        },
        "systematics_note": ("Shape statistical error is negligible (~%.2e). The "
                             "NORMALISATION carries the ENDF σ_nγ/σ_np systematic "
                             "(~25%%, overall scale) plus the α_ipc=%.1e capt→IPC and "
                             "BR_X17=%.1e assumptions — rescale via the scalars above."
                             % (1.0/np.sqrt(Nnp_ingate), a_ipc, br)),
        "usage": "Load with scripts/ipc_time_spectrum.py: density(t), cdf_at(t), "
                 "sample(n) -> unbinned arrival times, sample_energy(n).",
    }

    # ── save: npz (full), json (portable), two csv ────────────────────────────
    np.savez(OUT / "ipc_ingate_spectrum.npz",
             t_ms=tf, dNdt_ipc_per_pulse_per_ms=dNdt, cdf=cdf, E_eV=Ef,
             bin_E_lo=E_lo[ingate], bin_E_hi=np.minimum(E_hi, E_GATE)[ingate],
             bin_t_lo=t_of_E(np.minimum(E_hi, E_GATE)[ingate]),
             bin_t_hi=t_of_E(E_lo[ingate]),
             bin_ipc_per_pulse=ipc_bin[ingate],
             bin_np_counts=(h_np * frac_ig)[ingate],
             bin_ratio=rat_bin[ingate],
             ipc_per_pulse_ingate=ipc_ingate,
             ipc_per_day_ingate=ipc_ingate * PULSES_DAY,
             captures_per_pulse_ingate=cap_ingate,
             alpha_ipc=a_ipc, br_x17=br,
             pulses_per_day=PULSES_DAY, flight_ms=FLIGHT_MS, gate_ms=GATE_MS,
             meta_json=json.dumps(meta))

    bins_tbl = []
    for i in np.where(ingate)[0]:
        ehi = min(E_hi[i], E_GATE)
        bins_tbl.append({
            "E_lo_eV": E_lo[i], "E_hi_eV": ehi,
            "t_lo_ms": t_of_E(ehi), "t_hi_ms": t_of_E(E_lo[i]),
            "np_counts": (h_np[i] * frac_ig[i]),
            "ratio_ng_np": rat_bin[i],
            "ipc_per_pulse": ipc_bin[i],
            "ipc_per_day": ipc_bin[i] * PULSES_DAY,
            "ipc_per_pulse_stat_err": ipc_bin[i] / np.sqrt(max(h_np[i]*frac_ig[i], 1.0)),
        })
    meta_out = dict(meta)
    meta_out["bins"] = bins_tbl
    meta_out["fine_curve_note"] = ("full fine grid in the .npz/.csv; a 200-pt "
                                   "downsample is inlined here for portability.")
    ds = np.linspace(0, len(tf) - 1, 200).astype(int)
    meta_out["fine_downsampled"] = {"t_ms": tf[ds].tolist(),
                                    "dNdt_ipc_per_pulse_per_ms": dNdt[ds].tolist(),
                                    "cdf": cdf[ds].tolist()}
    (OUT / "ipc_ingate_spectrum.json").write_text(json.dumps(meta_out, indent=1))

    with open(OUT / "ipc_ingate_spectrum_fine.csv", "w") as f:
        f.write("# in-gate IPC production spectrum (reweighted). "
                "t_ms=arrival time; dNdt in IPC pairs/pulse/ms; "
                "cdf=fraction of in-gate IPC with arrival<=t; E_eV=neutron energy\n")
        f.write("# integral(dNdt dt)= %.6e IPC/pulse in-gate; ×%.4g pulses/day; "
                "×%.3g capt→IPC; ×%.3g X17/IPC\n" % (ipc_ingate, PULSES_DAY, a_ipc, br))
        f.write("t_ms,dNdt_ipc_per_pulse_per_ms,cdf,E_eV\n")
        for a, b, c, e in zip(tf, dNdt, cdf, Ef):
            f.write(f"{a:.6e},{b:.6e},{c:.6e},{e:.6e}\n")

    with open(OUT / "ipc_ingate_spectrum_bins.csv", "w") as f:
        f.write("# raw reweighted histogram (10 bins/decade). "
                "Fold this for exact per-bin content; errors are Poisson on (n,p)t.\n")
        f.write("E_lo_eV,E_hi_eV,t_lo_ms,t_hi_ms,np_counts,ratio_ng_np,"
                "ipc_per_pulse,ipc_per_day,ipc_per_pulse_stat_err\n")
        for r in bins_tbl:
            f.write(f"{r['E_lo_eV']:.6e},{r['E_hi_eV']:.6e},{r['t_lo_ms']:.6e},"
                    f"{r['t_hi_ms']:.6e},{r['np_counts']:.6e},{r['ratio_ng_np']:.6e},"
                    f"{r['ipc_per_pulse']:.6e},{r['ipc_per_day']:.6e},"
                    f"{r['ipc_per_pulse_stat_err']:.6e}\n")

    # ── the requested lin-lin plot (t > 1 ms only) ─────────────────────────────
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ax.fill_between(tf, dNdt, color="#1f77b4", alpha=0.18)
    ax.plot(tf, dNdt, color="#1f77b4", lw=2, label="reweighted (fine, unbinned-equiv.)")
    # raw per-bin points as a density (ipc/pulse ÷ Δt_bin) for validation
    bt_lo = t_of_E(np.minimum(E_hi, E_GATE)[ingate]); bt_hi = t_of_E(E_lo[ingate])
    bt_c  = 0.5 * (bt_lo + bt_hi); bdens = ipc_bin[ingate] / (bt_hi - bt_lo)
    m = ipc_bin[ingate] > 0
    ax.plot(bt_c[m], bdens[m], "o", ms=3.5, color="#d62728",
            label="raw bins (10/decade)")
    # two physical features: the epithermal edge at the gate, and the thermal peak
    it = np.argmax(np.where(tf > 3.0, dNdt, 0.0))         # thermal peak (t>3 ms)
    ax.axvline(tf[it], color="0.5", ls=":", lw=1)
    ax.annotate(f"thermal peak ≈ {tf[it]:.1f} ms\n(E ≈ {E_of_t(tf[it])*1e3:.0f} meV)",
                (tf[it], dNdt[it]), textcoords="offset points", xytext=(10, 6),
                fontsize=9, va="bottom")
    ax.annotate("epithermal edge\n(gate opens, E≈2 eV)", (tf[0], dNdt[0]),
                textcoords="offset points", xytext=(14, -4), fontsize=8.5, va="top",
                color="0.35")
    ax.set_xlim(GATE_MS, tf.max())
    ax.set_ylim(0, dNdt.max() * 1.12)
    ax.set_xlabel("neutron arrival time  t [ms]   (gate: t > 1 ms)")
    ax.set_ylabel("IPC pairs / pulse / ms")
    ax.set_title("In-gate IPC production spectrum (linear, reweighted)")
    ax.grid(True, alpha=0.25)
    sec = ax.secondary_yaxis("right",
        functions=(lambda y: y * PULSES_DAY, lambda y: y / PULSES_DAY))
    sec.set_ylabel("IPC pairs / day / ms")
    ax.legend(loc="upper right", fontsize=9,
              title=f"∫ = {ipc_ingate:.2e} IPC/pulse\n"
                    f"    = {ipc_ingate*PULSES_DAY:.2f} IPC/day")
    fig.tight_layout()
    fig.savefig(OUT / "ipc_vs_time_ingate_linlin.pdf", bbox_inches="tight")
    fig.savefig(OUT / "ipc_vs_time_ingate_linlin.png", dpi=140, bbox_inches="tight")

    print(f"in-gate IPC/pulse   : {ipc_ingate:.4e}  ({ipc_ingate*PULSES_DAY:.2f} IPC/day)")
    print(f"captures/pulse      : {cap_ingate:.4e}")
    print(f"eff (n,p)t counts   : {Nnp_ingate:.3e}  -> shape stat err {100/np.sqrt(Nnp_ingate):.3f}%")
    print(f"thermal peak time   : {tf[it]:.2f} ms   (E={E_of_t(tf[it])*1e3:.1f} meV)")
    print(f"time span (in-gate) : {tf.min():.2f} – {tf.max():.2f} ms")
    print(f"\nwrote -> {OUT}/")
    for fn in ["ipc_ingate_spectrum.npz", "ipc_ingate_spectrum.json",
               "ipc_ingate_spectrum_fine.csv", "ipc_ingate_spectrum_bins.csv",
               "ipc_vs_time_ingate_linlin.pdf"]:
        print(f"   {fn}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--npz", default=str(NPZ),
                    help="capture-scan npz (default: June sub-keV)")
    ap.add_argument("--outdir", default=str(OUT))
    a = ap.parse_args()
    main(npz_path=Path(a.npz), out_dir=Path(a.outdir))
