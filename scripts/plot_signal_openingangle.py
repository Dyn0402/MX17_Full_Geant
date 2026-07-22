#!/usr/bin/env python3
"""plot_signal_openingangle.py — first-look figures for the Al-pair background
analysis (Task 2 + Task 4a interim). Uses analysis/al_pair/signal_openingangle.npz
(X17/IPC truth from the pairs tree + low-stats Al conv pairs from the trajdump).

X17/IPC opening angles are truth from the generator. NOTE: the generator's IPC
is a 1/Mee mass spectrum + isotropic two-body decay + boost — it carries NO
multipole (E0/M1/E2) angular correlation, so the IPC curve here is a modelling
baseline, not a first-principles IPC angular distribution.

The Al curve is the LOW-STATS trajdump cross-check (n~22), and its angle is
biased UP by first-step multiple scattering — the clean high-stats truth comes
from the ConvPairTree campaign (neutrons_convpair_2cm).
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load("analysis/al_pair/signal_openingangle.npz")
x17, ipc, al = d["x17_theta"], d["ipc_theta"], d["al_theta"]

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))

# ── Panel A: opening-angle PDFs ─────────────────────────────────────────────
bins = np.linspace(0, 180, 91)
ax[0].hist(x17, bins=bins, density=True, histtype="stepfilled", alpha=0.35,
           color="C3", label=f"X17 (m=16.8 MeV)  med {np.median(x17):.0f}°")
ax[0].hist(ipc, bins=bins, density=True, histtype="step", lw=2.0,
           color="C0", label=f"IPC (1/Mee model)  med {np.median(ipc):.0f}°")
if len(al):
    ax[0].hist(al, bins=np.linspace(0, 180, 31), density=True, histtype="step",
               lw=1.6, color="0.35", ls="--",
               label=f"Al pair (trajdump n={len(al)}, MSC-biased)")
ax[0].set_xlabel("e⁺e⁻ opening angle  [deg]")
ax[0].set_ylabel("normalized / bin")
ax[0].set_title("Truth opening-angle distributions (infinite-stats shapes)")
ax[0].legend(fontsize=8.5, loc="upper left")
ax[0].set_xlim(0, 180)

# ── Panel B: individual lepton KE ───────────────────────────────────────────
x17_min = np.minimum(d["x17_em"], d["x17_ep"])
ipc_min = np.minimum(d["ipc_em"], d["ipc_ep"])
kb = np.linspace(0, 21, 106)
ax[1].hist(x17_min, bins=kb, density=True, histtype="step", lw=2, color="C3",
           label="X17: softer lepton KE")
ax[1].hist(ipc_min, bins=kb, density=True, histtype="step", lw=2, color="C0",
           label="IPC: softer lepton KE")
if len(d["al_emin_ke"]):
    ax[1].hist(d["al_emin_ke"], bins=np.linspace(0, 8, 33), density=True,
               histtype="step", lw=1.6, color="0.35", ls="--",
               label=f"Al pair: softer lepton KE (n={len(d['al_emin_ke'])})")
ax[1].axvspan(0, 3.0, color="0.85", alpha=0.5, zorder=0)
ax[1].set_xlabel("softer-lepton kinetic energy  [MeV]")
ax[1].set_ylabel("normalized / bin")
ax[1].set_title("Individual-lepton energy (a possible track-level handle)")
ax[1].legend(fontsize=8.5)
ax[1].set_xlim(0, 21)

fig.suptitle("Al pair-production background vs IPC / X17 — opening angle & lepton energy",
             fontsize=12, y=1.02)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"analysis/al_pair/signal_openingangle.{ext}", bbox_inches="tight", dpi=140)
print("wrote analysis/al_pair/signal_openingangle.png/.pdf")
print(f"X17 opening: med={np.median(x17):.1f}  IQR=[{np.percentile(x17,25):.0f},{np.percentile(x17,75):.0f}]")
print(f"IPC opening: med={np.median(ipc):.1f}  IQR=[{np.percentile(ipc,25):.0f},{np.percentile(ipc,75):.0f}]")
print(f"IPC fraction above 100 deg (X17 region): {(ipc>100).mean()*100:.2f}%")
print(f"X17 fraction below 40 deg (IPC region):  {(x17<40).mean()*100:.2f}%")
