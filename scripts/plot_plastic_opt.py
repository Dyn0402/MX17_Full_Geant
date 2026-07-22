#!/usr/bin/env python3
"""plot_plastic_opt.py — optimize the PLASTIC trigger threshold to select real MM
tracks (IPC) and reject double-Compton (photon-only, no-track) triggers, given
the DAQ can read only ~4 events on the thermal peak.

Inputs (from digest_plastic_opt.py): analysis/al_pair/plastic_opt_{signal,thermal}.npz
Per (event,arm): S,P (SiPM/plastic max-bar MeV), MM (DriftGas+AmpGas charged MeV,
the TRACK label), EP (e+ in scint), typ (0=X17,1=IPC), isAl.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MIP_S, MIP_P = 0.458, 4.334 * 2.0 / 2.5      # SiPM, plastic(2cm scaled) MeV/MIP
S_THR   = 0.5 * MIP_S                          # fixed SiPM leg threshold [MeV]
MM_THR  = 0.005                                # MM charged edep for a "track" [MeV]
N_PULSE = 4.284e6                              # thermal n/pulse (in-gate window)
# DAQ readout budget = events readable per spill in the thermal band.
# 2026-07-22: 10 GbE upgrade + readout start pushed to 1 ms → ~24 ev/spill
# (docs/network_upgrade_10g/04_bandwidth_model.md, IPD 10; was ~4 pre-upgrade).
DAQ_BUDGET = 24.0                              # readable events / thermal pulse

sig = np.load("analysis/al_pair/plastic_opt_signal.npz")
th  = np.load("analysis/al_pair/plastic_opt_thermal.npz")

def per_event(d):
    return d["ev"], d["arm"], d["S"], d["P"], d["MM"], d["EP"], d["typ"], d["isAl"], int(d["n_events"])

sev, sarm, sS, sP, sMM, sEP, styp, _, s_nev = per_event(sig)
tev, tarm, tS, tP, tMM, tEP, _, tAl, t_nev = per_event(th)
w_th = N_PULSE / t_nev                          # thermal count -> /pulse
n_x17 = s_nev / 2.0; n_ipc = s_nev / 2.0        # 50/50 generated

def ev_any(ev, mask):
    return np.unique(ev[mask])

# ── plastic-edep spectra of leg arms, by class (SiPM leg required) ───────────
sig_leg = sS >= S_THR
th_leg  = tS >= S_THR
ipc_track = sig_leg & (styp == 1) & (sMM >= MM_THR)     # IPC real-track legs
x17_track = sig_leg & (styp == 0) & (sMM >= MM_THR)
bg_dc = th_leg & tAl & (tMM < MM_THR) & (tEP <= 0)      # Al double-Compton (no track)
bg_tk = th_leg & tAl & (tMM >= MM_THR)                  # Al track legs (pair-like)
print(f"MM_THR={MM_THR*1e3:.0f} keV | SiPM leg S>= {S_THR*1e3:.0f} keV")
print(f"leg arms: IPC-track={ipc_track.sum()}  Al double-Compton={bg_dc.sum()}  Al track={bg_tk.sum()}")
print(f"plastic edep [MeV] medians: IPC-track={np.median(sP[ipc_track]):.2f}  "
      f"Al-DC={np.median(tP[bg_dc]):.2f}  Al-track={np.median(tP[bg_tk]):.2f}")

# ── threshold scan (SiPM fixed at 0.5 MIP; scan plastic b) ───────────────────
b_grid = np.linspace(0.05, 8.0, 160)           # plastic threshold [MeV]
ipc_eff1 = []; ipc_eff2 = []; x17_eff1 = []
bg_all = []; bg_notrack = []; bg_track = []
for b in b_grid:
    sl = sig_leg & (sP >= b)                    # signal leg arms
    tl = th_leg & (tP >= b)                     # thermal leg arms
    # signal per-event efficiency
    def eff(typ_code, ngen, need=1):
        m = sl & (styp == typ_code)
        ev = sev[m]
        if need == 1:
            return len(np.unique(ev)) / ngen
        u, cnt = np.unique(ev, return_counts=True)
        return (cnt >= 2).sum() / ngen
    ipc_eff1.append(eff(1, n_ipc, 1)); ipc_eff2.append(eff(1, n_ipc, 2))
    x17_eff1.append(eff(0, n_x17, 1))
    # background per-event trigger load (any leg), and track vs no-track
    ev_trig = np.unique(tev[tl])
    ev_track = np.unique(tev[tl & (tMM >= MM_THR)])   # event has >=1 track-leg
    n_trig = len(ev_trig); n_track = len(ev_track)
    bg_all.append(n_trig * w_th)
    bg_track.append(n_track * w_th)
    bg_notrack.append((n_trig - n_track) * w_th)
ipc_eff1 = np.array(ipc_eff1); ipc_eff2 = np.array(ipc_eff2); x17_eff1 = np.array(x17_eff1)
bg_all = np.array(bg_all); bg_notrack = np.array(bg_notrack); bg_track = np.array(bg_track)

# Operationally-correct FOM: IPC is ~1e-4/pulse — far rarer than the ~200/pulse
# background trigger stream. When the DAQ saturates (bg>budget) it reads a random
# 4, so a triggered IPC is read with prob min(1, budget/bg). Maximize the IPC
# read-yield = eff × min(1, budget/bg).
read_frac = np.minimum(1.0, DAQ_BUDGET / bg_all)
ipc_read = ipc_eff1 * read_frac                    # ∝ IPC events read/pulse
cur = np.argmin(np.abs(b_grid - 0.5 * MIP_P))      # current 0.5 MIP plastic
best = int(np.argmax(ipc_read))
b_budget = b_grid[np.argmin(np.abs(bg_all - DAQ_BUDGET))]

print(f"\n current 0.5 MIP plastic = {0.5*MIP_P:.2f} MeV:")
print(f"   IPC eff(>=1)={ipc_eff1[cur]*100:.1f}%  bg trig/pulse={bg_all[cur]:.0f}  "
      f"(double-Compton {bg_notrack[cur]:.0f} + track {bg_track[cur]:.0f})  "
      f"DAQ read-frac={read_frac[cur]*100:.1f}%  -> IPC read-yield (rel)=1.00")
print(f"\n RECOMMEND plastic ~{b_grid[best]:.2f} MeV = {b_grid[best]/MIP_P:.2f} MIP "
      f"(max IPC read-yield):")
print(f"   IPC eff(>=1)={ipc_eff1[best]*100:.1f}%  bg trig/pulse={bg_all[best]:.1f}  "
      f"(double-Compton {bg_notrack[best]:.1f} + track {bg_track[best]:.1f})")
print(f"   DAQ read-frac={read_frac[best]*100:.0f}%  -> IPC read-yield ×{ipc_read[best]/ipc_read[cur]:.1f} vs current")
print(f"\n background reaches the DAQ budget ({DAQ_BUDGET:.0f}/pulse) at plastic "
      f"~{b_budget:.2f} MeV = {b_budget/MIP_P:.2f} MIP")
ir = best; brec = b_grid[best]

# ── sensitivity to the DAQ budget (the number that just changed) ─────────────
print(f"\n== optimum plastic threshold vs DAQ budget (bg rate 0.5 MIP = {bg_all[cur]:.0f}/pulse) ==")
print(f"{'budget/pulse':>12} {'opt thr MeV':>11} {'opt MIP':>8} {'IPC eff':>8} {'IPC read ×vs 0.5MIP@thatbudget':>32}")
for bud in (4, 12, 24, 48, 96):
    rf = np.minimum(1.0, bud / bg_all)
    read = ipc_eff1 * rf
    bi = int(np.argmax(read))
    base = ipc_eff1[cur] * min(1.0, bud / bg_all[cur])
    print(f"{bud:>12d} {b_grid[bi]:>11.2f} {b_grid[bi]/MIP_P:>8.2f} {ipc_eff1[bi]*100:>7.1f}% "
          f"{read[bi]/base:>30.1f}×")

# ── figure ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 3, figsize=(18, 4.9))
bins = np.linspace(0, 12, 61)
for m, P, col, lab in [(ipc_track, sP, "C3", "IPC real-track leg"),
                       (bg_tk, tP, "C0", "Al track leg (pair)"),
                       (bg_dc, tP, "0.3", "Al double-Compton leg (no track)")]:
    ax[0].hist(P[m], bins=bins, density=True, histtype="step", lw=2, color=col, label=lab)
ax[0].axvline(0.5*MIP_P, color="k", ls=":", lw=1, label=f"current 0.5 MIP={0.5*MIP_P:.1f} MeV")
ax[0].axvline(brec, color="green", ls="--", lw=1.5, label=f"recommend {brec:.1f} MeV")
ax[0].set(xlabel="plastic max-bar edep [MeV]", ylabel="norm / bin",
          title="(A) Plastic energy by class")
ax[0].legend(fontsize=8)

axb = ax[1]; axr = axb.twinx()
axb.plot(b_grid, ipc_eff1*100, "C3-", lw=2, label="IPC eff (≥1 leg)")
axb.plot(b_grid, ipc_eff2*100, "C3--", lw=1.6, label="IPC eff (≥2 legs)")
axr.plot(b_grid, bg_notrack, "-", color="0.3", lw=2, label="double-Compton trig/pulse")
axr.plot(b_grid, bg_track, "-", color="C0", lw=1.6, label="track bg trig/pulse")
axr.axhline(DAQ_BUDGET, color="purple", ls=":", lw=1.2)
axr.set_yscale("log")
axb.axvline(0.5*MIP_P, color="k", ls=":", lw=1); axb.axvline(brec, color="green", ls="--", lw=1.5)
axb.set(xlabel="plastic threshold [MeV]", ylabel="signal efficiency [%]",
        title="(B) Efficiency & background vs threshold")
axr.set_ylabel("background triggers / pulse")
axb.legend(fontsize=8, loc="lower left"); axr.legend(fontsize=8, loc="upper right")

ax[2].plot(b_grid, ipc_read / ipc_read[cur], "-", color="0.15", lw=2.2,
           label="IPC read-yield (rel. to 0.5 MIP)")
ax[2].axvline(0.5*MIP_P, color="k", ls=":", lw=1, label=f"current 0.5 MIP")
ax[2].axvline(brec, color="green", ls="--", lw=1.5, label=f"optimum {brec:.1f} MeV = {brec/MIP_P:.1f} MIP")
ax[2].axvline(b_budget, color="purple", ls=":", lw=1.2, label=f"bg=DAQ budget @ {b_budget:.1f} MeV")
ax[2].set(xlabel="plastic threshold [MeV]", ylabel="IPC events read / pulse (relative)",
          title="(C) DAQ throughput: IPC read-yield vs threshold")
ax[2].legend(fontsize=8)
fig.suptitle(f"Plastic threshold optimization — de-saturate the DAQ: "
             f"raise 0.5→{brec/MIP_P:.1f} MIP (~{brec:.1f} MeV) records ×{ipc_read[best]/ipc_read[cur]:.0f} more IPC",
             y=1.02, fontsize=12)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"analysis/al_pair/plastic_opt.{ext}", bbox_inches="tight", dpi=140)
print("\nwrote analysis/al_pair/plastic_opt.png/.pdf")
