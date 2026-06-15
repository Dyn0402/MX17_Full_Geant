#!/usr/bin/env python3
"""
make_e0_final_metric.py

The "final-metric" deliverable for the E0 pair-channel thread: X17 and IPC
(e+e-) pairs per pulse / per day as a function of neutron energy, with the
gamma-dark E0 channel included, folded through detector acceptance AND the
trigger/readout time window.

It combines four ingredients (see docs/e0_branch/):

  1. FORMATION x DECAY shapes (Part 1 x Part 2 of the E0 note):
       term 1 (M1/E1): N1(E) = rad_per_pulse(E) * alpha_IPC      -> peaks at MeV
       term 2 (E0):    N0(E) = N_np(E) * f * ratio_thermal * g0+(E) -> sub-keV
     term 1 uses the DIRECT 3He(n,g)4He counts per decade from the 5e8-neutron
     campaign (analysis/mev/mev_rates.json) so it reproduces the committed
     mev_note headline (32.7 X17/day produced, alpha=2.1e-3); term 2 reuses the
     sub-threshold-0+ construction of make_e0_pair_yield_fig.py.

  2. ABSOLUTE NORMALISATION: pulses/day, alpha_IPC, BR_X17 (mev_rates.json /
     mev_note); E0 strength f = sigma_E0/sigma_M1 carried as a band 1e-3..1e-2.

  3. DETECTOR ACCEPTANCE: flat MM-double geometric acceptance from the 10M-event
     pairs_v2 at-rest study (docs/angular_resolution): 19.6% (X17), 23.6% (IPC).

  4. TRIGGER / READOUT TIME WINDOW (the part that decides whether E0 is even
     seen).  At L = 19.5 m (EAR2), TOF(E_n) = L / v:
       - FLASH readout, ~10 us/pulse  <=>  E_n >~ 20 keV  -> catches M1/E1,
         MISSES the sub-keV E0.
       - + THERMAL trigger             ->  recovers the sub-keV decades where
         E0 lives -> this is the scenario that makes E0 (and E0->X17) visible.

E0 -> X17:  a 0+ -> 0+ transition is gamma-dark only because the photon is
massless (transverse only).  X17 is massive, so a VECTOR (1-) or SCALAR (0+)
X17 CAN be emitted in the monopole transition (pseudoscalar/axial cannot).  We
therefore show an E0->X17 signal band, with its branching carried as a separate
benchmark (eta_E0, default = the same 2.5% X17/pair as the M1 line, flagged).

Outputs:
  docs/report/figs/fig_e0_final_metric.{pdf,png}   (production vs recorded)
  analysis/e0/final_metric.json                    (the money table numbers)
"""
from pathlib import Path
import json
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RATES = ROOT / "analysis" / "mev" / "mev_rates.json"
H5 = ROOT / "data" / "He3.h5"
FIGOUT = ROOT / "docs" / "report" / "figs"
JSONOUT = ROOT / "analysis" / "e0" / "final_metric.json"

# --- normalisations -------------------------------------------------------
ALPHA_ANCHOR = 3.5e-3        # high-energy IPC coeff (12C/8Be-anchored)
ALPHA_TABLE = 2.1e-3         # Alberto's rate-table value (mev_note headline)
BR_X17 = 0.025               # X17 per IPC pair (M1 line); benchmark for E0 too
F_E0 = [1e-3, 3e-3, 1e-2]    # sigma_E0/sigma_M1 bracket (lo, mid, hi)
ACC_X17 = 0.196              # MM-double geometric acceptance, X17
ACC_IPC = 0.236              # MM-double geometric acceptance, IPC
# trigger / TOF
L_FLIGHT = 19.5              # m, EAR2
MN_MEV = 939.565
C_MS = 2.998e8
T_FLASH_US = 10.0            # flash-readout window length
# sub-threshold 0+ resonance (TUNL): E_x=20.21, threshold 20.578 -> E_r=-0.368
ER_0P, G_0P = 20.21 - 20.578, 0.50
TEMP = "294K"


def tof_us(E_eV):
    """neutron TOF over L_FLIGHT, in us (non-relativistic, fine below ~10 MeV)."""
    E = np.asarray(E_eV, float) * 1e-6  # MeV
    beta = np.sqrt(2 * E / MN_MEV)
    return L_FLIGHT / (beta * C_MS) * 1e6


def E_at_tof(t_us):
    """neutron energy [eV] whose TOF over L_FLIGHT equals t_us."""
    beta = L_FLIGHT / (t_us * 1e-6 * C_MS)
    return 0.5 * MN_MEV * beta**2 * 1e6


def endf_ratio(E_MeV):
    """sigma_ng/sigma_np vs energy from ENDF/B-VIII.0 (data/He3.h5)."""
    with h5py.File(H5, "r") as f:
        g = f["He3"]
        def xs(rx):
            E = g[f"energy/{TEMP}"][...]
            ds = g[f"reactions/{rx}/{TEMP}/xs"]
            i0 = ds.attrs.get("threshold_idx", 0)
            return E[i0:i0 + ds.shape[0]] * 1e-6, ds[...]
        Eg, ng = xs("reaction_102")
        Ep, npx = xs("reaction_103")
    def il(x, xp, fp):
        return 10 ** np.interp(np.log10(x), np.log10(xp),
                               np.log10(np.clip(fp, 1e-40, None)))
    return il(E_MeV, Eg, ng) / il(E_MeV, Ep, npx)


def g0plus(E_cm_MeV):
    num = ER_0P**2 + (G_0P / 2)**2
    return num / ((E_cm_MeV - ER_0P)**2 + (G_0P / 2)**2)


def flash_fraction(Elo, Ehi):
    """fraction of a decade [Elo,Ehi] (eV) that arrives within the 10us flash
    window, i.e. with E_n >= E(10us) ~ 20 keV.  Log-uniform within the decade."""
    Ecut = E_at_tof(T_FLASH_US)
    if Ehi <= Ecut:
        return 0.0
    if Elo >= Ecut:
        return 1.0
    return (np.log10(Ehi) - np.log10(Ecut)) / (np.log10(Ehi) - np.log10(Elo))


def main():
    d = json.load(open(RATES))
    ppd = d["pulses_per_day"]
    rows = d["decades"]
    Elo = np.array([r["E_lo_eV"] for r in rows])
    Ehi = np.array([r["E_hi_eV"] for r in rows])
    Ecen = np.sqrt(Elo * Ehi)
    Ecen_MeV = Ecen * 1e-6
    rad_pp = np.array([r["rad_per_pulse"] for r in rows])              # direct (n,g)/pulse
    Np = np.array([r["beam_per_pulse"] * r["np_frac_of_beam"] for r in rows])

    ratio = endf_ratio(Ecen_MeV)
    ratio_th = endf_ratio(np.array([2.53e-8]))[0]
    Ecm_MeV = 0.75 * Ecen_MeV

    # ---- production, pairs / pulse / decade --------------------------------
    term1 = rad_pp * ALPHA_ANCHOR                                      # M1/E1 IPC pairs
    term1_tab = rad_pp * ALPHA_TABLE
    e0 = {f: Np * f * ratio_th * g0plus(Ecm_MeV) for f in F_E0}        # E0 pairs (=transitions)

    # ---- trigger fractions per decade --------------------------------------
    fflash = np.array([flash_fraction(lo, hi) for lo, hi in zip(Elo, Ehi)])
    fall = np.ones_like(fflash)                                        # flash + thermal trigger

    Ecut = E_at_tof(T_FLASH_US)

    def totals(per_pulse_decade, acc, tofmask):
        """-> per-day total over decades, given acceptance and a TOF fraction."""
        return float((per_pulse_decade * tofmask).sum() * acc * ppd)

    # IPC background pairs/day: M1/E1 + E0
    res = {
        "meta": {
            "alpha_anchor": ALPHA_ANCHOR, "alpha_table": ALPHA_TABLE,
            "BR_X17": BR_X17, "acc_X17": ACC_X17, "acc_IPC": ACC_IPC,
            "L_flight_m": L_FLIGHT, "t_flash_us": T_FLASH_US,
            "E_cut_flash_eV": float(Ecut), "f_E0": F_E0,
            "pulses_per_day": ppd, "ratio_thermal_ng_np": float(ratio_th),
            "note": "term1 = direct G4 (n,g) counts x alpha; term2 (E0) = "
                    "N_np x f x ratio_thermal x g0+(E). recorded = produced x "
                    "acceptance x TOF-window fraction.",
        },
        "produced_per_day": {}, "recorded_flash": {}, "recorded_thermal": {},
    }

    # --- IPC pairs (background) ---
    for tag, mask, store in [("flash", fflash, "recorded_flash"),
                             ("thermal", fall, "recorded_thermal")]:
        acc = ACC_IPC
        res[store]["ipc_M1E1"] = totals(term1, acc, mask)
        res[store]["ipc_E0"] = {f"f={f:g}": totals(e0[f], acc, mask) for f in F_E0}
    res["produced_per_day"]["ipc_M1E1"] = float(term1.sum() * ppd)
    res["produced_per_day"]["ipc_M1E1_table_alpha"] = float(term1_tab.sum() * ppd)
    res["produced_per_day"]["ipc_E0"] = {f"f={f:g}": float(e0[f].sum() * ppd) for f in F_E0}

    # --- X17 signal:  M1 line (rad x alpha x BR) + E0->X17 (E0 x eta_E0) ---
    x17_M1 = term1 * BR_X17
    x17_E0 = {f: e0[f] * BR_X17 for f in F_E0}     # eta_E0 = BR_X17 benchmark
    res["produced_per_day"]["x17_M1"] = float(x17_M1.sum() * ppd)
    res["produced_per_day"]["x17_M1_table_alpha"] = float((term1_tab * BR_X17).sum() * ppd)
    res["produced_per_day"]["x17_E0"] = {f"f={f:g}": float(x17_E0[f].sum() * ppd) for f in F_E0}
    for tag, mask, store in [("flash", fflash, "recorded_flash"),
                             ("thermal", fall, "recorded_thermal")]:
        res[store]["x17_M1"] = totals(x17_M1, ACC_X17, mask)
        res[store]["x17_E0"] = {f"f={f:g}": totals(x17_E0[f], ACC_X17, mask) for f in F_E0}

    JSONOUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(JSONOUT, "w"), indent=1)

    # ======================= figure ========================================
    def nz(a):
        """mask non-positive entries to NaN so log plots don't draw drop-lines."""
        a = np.asarray(a, float).copy()
        a[~(a > 0)] = np.nan
        return a

    fig, (axB, axS) = plt.subplots(1, 2, figsize=(13.5, 5.8), sharex=True, sharey=True)

    for ax, title in [(axB, "IPC (e$^+$e$^-$) background pairs"),
                      (axS, "X17 signal")]:
        # trigger-window shading
        ax.axvspan(1e-3, Ecut, color="C0", alpha=0.08, zorder=0)
        ax.axvspan(Ecut, 1e8, color="gold", alpha=0.12, zorder=0)
        ax.axvline(Ecut, color="0.5", ls=":", lw=1)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"neutron energy  $E_n$  [eV]")
        ax.set_xlim(1e-2, 1e8); ax.set_ylim(1e-3, 1e3)
        ax.grid(alpha=0.3, which="both")
        ax.set_title(title)

    # text labels for the two trigger regions (once, on the left panel)
    axB.text(2e-1, 4e2, "thermal trigger\nrecovers this", ha="center",
             fontsize=8.5, color="C0")
    axB.text(2e6, 4e2, "flash 10 µs\nreadout", ha="center", fontsize=8.5,
             color="0.45")

    # --- background panel ---
    axB.plot(Ecen, nz(term1 * ppd), "o-", color="C0", ms=6, lw=2.2,
             label=r"M1/E1 IPC produced ($\alpha=3.5\times10^{-3}$)")
    axB.fill_between(Ecen, nz(e0[F_E0[0]] * ppd), nz(e0[F_E0[2]] * ppd),
                     color="C3", alpha=0.22,
                     label=r"E0 produced, $f=10^{-3}$–$10^{-2}$")
    axB.plot(Ecen, nz(e0[F_E0[1]] * ppd), "--", color="C3", lw=1.8,
             label=r"E0 produced, $f=3\times10^{-3}$")
    # recorded (flash gate) overlay -- only M1/E1 survives; E0 -> ~0 (annotated)
    axB.plot(Ecen, nz(term1 * ACC_IPC * fflash * ppd), "s:", color="C0", ms=4,
             lw=1.3, alpha=0.85, label="M1/E1 recorded (flash×acc)")
    axB.annotate("E0 recorded in flash window ≈ 0\n(E0 is sub-keV; needs a thermal trigger)",
                 xy=(3e0, 1.3e-2), fontsize=7.8, color="C3", ha="center")
    axB.legend(fontsize=7.6, loc="lower right", ncol=1)

    # --- signal panel ---
    axS.plot(Ecen, nz(x17_M1 * ppd), "o-", color="C2", ms=6, lw=2.2,
             label=r"X17 from M1 produced")
    axS.fill_between(Ecen, nz(x17_E0[F_E0[0]] * ppd), nz(x17_E0[F_E0[2]] * ppd),
                     color="C1", alpha=0.25,
                     label="X17 from E0 produced\n"
                           r"($f=10^{-3}$–$10^{-2}$, vector/scalar X17)")
    axS.plot(Ecen, nz(x17_M1 * ACC_X17 * fflash * ppd), "s:", color="C2", ms=4,
             lw=1.3, alpha=0.85, label="X17(M1) recorded (flash×acc)")
    axS.legend(fontsize=7.8, loc="lower right")

    axB.set_ylabel("pairs (or X17) produced / recorded per day / decade")
    fig.suptitle("Final metric: X17 / IPC pairs vs neutron energy, with the "
                 "γ-dark E0 channel and the trigger window\n"
                 f"flash 10 µs readout reaches $E_n\\gtrsim${Ecut/1e3:.0f} keV; "
                 "E0 lives sub-keV → only a thermal trigger records it",
                 fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    FIGOUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIGOUT / f"fig_e0_final_metric.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

    # ======================= console summary ===============================
    def fmt(x): return f"{x:8.2f}"
    print(f"E_cut(flash 10us) = {Ecut/1e3:.1f} keV   ratio_thermal = {ratio_th:.2e}")
    print(f"pulses/day = {ppd:.0f}   acc X17/IPC = {ACC_X17}/{ACC_IPC}\n")
    print("                         produced/day   rec.flash/day   rec.+thermal/day")
    print(f"IPC M1/E1 (a=3.5e-3) {fmt(res['produced_per_day']['ipc_M1E1'])}"
          f"   {fmt(res['recorded_flash']['ipc_M1E1'])}   {fmt(res['recorded_thermal']['ipc_M1E1'])}")
    print(f"IPC M1/E1 (a=2.1e-3) {fmt(res['produced_per_day']['ipc_M1E1_table_alpha'])}")
    for f in F_E0:
        k = f"f={f:g}"
        print(f"IPC E0  {k:11s} {fmt(res['produced_per_day']['ipc_E0'][k])}"
              f"   {fmt(res['recorded_flash']['ipc_E0'][k])}   {fmt(res['recorded_thermal']['ipc_E0'][k])}")
    print(f"X17 M1   (a=3.5e-3) {fmt(res['produced_per_day']['x17_M1'])}"
          f"   {fmt(res['recorded_flash']['x17_M1'])}   {fmt(res['recorded_thermal']['x17_M1'])}")
    print(f"X17 M1   (a=2.1e-3) {fmt(res['produced_per_day']['x17_M1_table_alpha'])}")
    for f in F_E0:
        k = f"f={f:g}"
        print(f"X17 E0  {k:11s} {fmt(res['produced_per_day']['x17_E0'][k])}"
              f"   {fmt(res['recorded_flash']['x17_E0'][k])}   {fmt(res['recorded_thermal']['x17_E0'][k])}")
    print(f"\nSaved -> {FIGOUT}/fig_e0_final_metric.[pdf,png]")
    print(f"Saved -> {JSONOUT}")


if __name__ == "__main__":
    main()
