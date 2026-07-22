#!/usr/bin/env python3
"""plot_convpair_truth.py — truth-level Al pair background vs IPC/X17, with REAL
high-stats Al conv pairs from the ConvPairTree campaign (analysis/al_pair/
convpair_truth.npz) overlaid on the signal (signal_openingangle.npz).

Panels: (A) normalized opening angle, (B) softer-lepton KE, (C) rate-weighted
opening angle per day (Al produced rate vs IPC/X17 in-gate) — the raw dominance.
"""
import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

s = np.load("analysis/al_pair/signal_openingangle.npz")
c = np.load("analysis/al_pair/convpair_truth.npz")
meta = json.load(open("analysis/al_pair/convpair_truth.json"))
cv = np.array([x.decode() if isinstance(x, bytes) else x for x in c["convvol"]])
al = cv == "He3Cap_Al"
al_theta = c["theta"][al]
al_emin  = c["emin"][al]
x17, ipc = s["x17_theta"], s["ipc_theta"]
X17_THRESH = 108.0

# in-gate signal rates (handoff / reweight_2cm); Al = produced rate (Task-3 MM
# acceptance handled separately in the reco analysis)
R_AL_DAY  = meta["al_pairs_per_day"]      # 5.95e6/day produced
R_IPC_DAY = 1.39
R_X17_DAY = 0.035

fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.7))

# ── A: normalized opening angle ─────────────────────────────────────────────
b = np.linspace(0, 180, 91)
ax[0].hist(x17, bins=b, density=True, histtype="stepfilled", alpha=0.35, color="C3",
           label=f"X17  med {np.median(x17):.0f}°")
ax[0].hist(ipc, bins=b, density=True, histtype="step", lw=2, color="C0",
           label=f"IPC (1/Mee model)  med {np.median(ipc):.0f}°")
ax[0].hist(al_theta, bins=b, density=True, histtype="step", lw=2, color="0.25",
           label=f"Al pair TRUTH (n={al.sum()})  med {np.median(al_theta):.0f}°")
ax[0].axvline(X17_THRESH, color="C3", ls=":", lw=1.2)
ax[0].set(xlabel="e⁺e⁻ opening angle [deg]", ylabel="normalized / bin", xlim=(0, 180),
          title="(A) Truth opening-angle shapes")
ax[0].legend(fontsize=8.5, loc="upper center")

# ── B: softer-lepton KE ─────────────────────────────────────────────────────
kb = np.linspace(0, 21, 106)
x17m = np.minimum(s["x17_em"], s["x17_ep"]); ipcm = np.minimum(s["ipc_em"], s["ipc_ep"])
ax[1].hist(x17m, bins=kb, density=True, histtype="step", lw=2, color="C3", label="X17")
ax[1].hist(ipcm, bins=kb, density=True, histtype="step", lw=2, color="C0", label="IPC")
ax[1].hist(al_emin, bins=kb, density=True, histtype="step", lw=2, color="0.25",
           label=f"Al pair TRUTH")
ax[1].axvspan(0, 4.0, color="0.85", alpha=0.6, zorder=0)
ax[1].text(2.0, ax[1].get_ylim()[1]*0.9, "Al lives here\n(X17 never does)",
           ha="center", va="top", fontsize=8, color="0.3")
ax[1].set(xlabel="softer-lepton kinetic energy [MeV]", ylabel="normalized / bin",
          xlim=(0, 21), title="(B) Softer-lepton energy (track-momentum handle)")
ax[1].legend(fontsize=9)

# ── C: rate-weighted opening angle per day ──────────────────────────────────
bw = b[1] - b[0]
for data, R, col, lab in [(al_theta, R_AL_DAY, "0.25", "Al pair (produced)"),
                          (ipc, R_IPC_DAY, "C0", "IPC (in-gate)"),
                          (x17, R_X17_DAY, "C3", "X17 (in-gate)")]:
    h, _ = np.histogram(data, bins=b, density=True)
    ax[2].step(b[:-1], h * R * bw, where="post", color=col, lw=2, label=f"{lab}  {R:.3g}/day")
ax[2].axvline(X17_THRESH, color="C3", ls=":", lw=1.2)
ax[2].set_yscale("log")
ax[2].set(xlabel="e⁺e⁻ opening angle [deg]", ylabel="pairs / day / bin", xlim=(0, 180),
          title="(C) Rate-weighted (produced Al vs in-gate signal)")
ax[2].legend(fontsize=8.5, loc="upper right")
ax[2].set_ylim(1e-4, 1e6)

fig.suptitle("Al pair-production background vs IPC / X17 — truth kinematics "
             f"(Al: {meta['al_pairs_per_day']:.2g}/day produced)", y=1.02, fontsize=12)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"analysis/al_pair/convpair_truth.{ext}", bbox_inches="tight", dpi=140)
print("wrote analysis/al_pair/convpair_truth.png/.pdf")
print(f"Al truth: median={np.median(al_theta):.1f}°  >108°={ (al_theta>108).mean()*100:.2f}%  "
      f"-> {(al_theta>108).mean()*R_AL_DAY:.3g}/day above X17 threshold (produced)")
print(f"Al softer-lepton >4 MeV: {(al_emin>4).mean()*100:.4f}%")
