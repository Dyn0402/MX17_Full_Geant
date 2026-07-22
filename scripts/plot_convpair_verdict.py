#!/usr/bin/env python3
"""plot_convpair_verdict.py — Task 5: the feasibility verdict. Combines the
detector-level (MM-reconstructed, MSC-smeared) opening angle of the Al pair
background (convpair_reco.npz) and of the signal (signal_reco.npz), normalizes
to pairs/day, and scans opening-angle (± soft-lepton momentum) selections to get
the residual Al under the X17 peak and the achievable S/B.

Rates (in-gate PRODUCTION, thermal gate):
  Al capsule conv : 5.95e6/day  (analyze_convpair_truth: He3Cap_Al produced)
  IPC             : 1.39 /day   (reweight_2cm)
  X17             : 0.035/day   (2.5% of IPC)
MM 2-track acceptance for each is measured from its own *_reco.npz.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R_AL, R_IPC, R_X17 = 5.95e6, 1.39, 0.035     # in-gate production /day

a = np.load("analysis/al_pair/convpair_reco.npz")
s = np.load("analysis/al_pair/signal_reco.npz")

# Al MM 2-track pairs
al_pair = a["n_mm"] == 2
al_acc  = al_pair.mean()
al_tr   = a["theta_r"][al_pair]
al_emin = a["emin_t"][al_pair]
# signal MM 2-track pairs
def sig(code):
    m = s["type"] == code
    p = m & (s["n_mm"] == 2)
    return p[m].mean(), s["theta_reco"][p], s["emin_truth"][p]
x17_acc, x17_tr, x17_emin = sig(0)
ipc_acc, ipc_tr, ipc_emin = sig(1)

# MM-pair rates /day = production × acceptance
RA, RI, RX = R_AL*al_acc, R_IPC*ipc_acc, R_X17*x17_acc
wA, wI, wX = RA/len(al_tr), RI/len(ipc_tr), RX/len(x17_tr)   # per-entry weight

print(f"MM 2-track acceptance:  Al={al_acc*100:.1f}%  IPC={ipc_acc*100:.1f}%  X17={x17_acc*100:.1f}%")
print(f"MM-pair rates /day:     Al={RA:.3g}  IPC={RI:.3g}  X17={RX:.3g}")

def survive(cut_lo, cut_hi, emin_cut):
    """rate/day of each species passing [cut_lo,cut_hi] deg AND both leptons >
    emin_cut MeV (soft-lepton handle; emin_cut=0 disables)."""
    def r(tr, emin, w):
        m = (tr >= cut_lo) & (tr <= cut_hi) & (emin >= emin_cut)
        return m.sum() * w
    return r(al_tr, al_emin, wA), r(ipc_tr, ipc_emin, wI), r(x17_tr, x17_emin, wX)

# X17 peak window (reco). Find the window that maximizes S/sqrt(B) roughly:
print("\n== residual under an X17 opening-angle window ==")
print(f"{'window':>12} {'soft-lep':>9} {'Al/day':>10} {'IPC/day':>9} {'X17/day':>9} {'Al:X17':>10}")
for lo, hi in [(108, 180), (110, 140), (100, 150)]:
    for ec in (0.0, 4.0):
        al_r, ipc_r, x17_r = survive(lo, hi, ec)
        ratio = al_r / x17_r if x17_r else np.inf
        print(f"{f'{lo}-{hi}°':>12} {f'>{ec:.0f}MeV' if ec else 'none':>9} "
              f"{al_r:>10.3g} {ipc_r:>9.3g} {x17_r:>9.3g} {ratio:>10.3g}")

# ── momentum-resolution scan (the crux) ─────────────────────────────────────
# TWO tracking-only energy handles (no calorimeter), each needs MM momentum:
#   (1) softer-lepton > 4 MeV  — Al softer is kinematically <3.35 MeV (7.72 line),
#       X17 softer >4 MeV: a thin ~0.65 MeV/18% gap → needs excellent resolution.
#   (2) TOTAL pair KE > 13 MeV — Al(7.72) total KE ~6.7 MeV; BOTH X17 and IPC
#       carry the full 20.58 transition → total KE ~19.6 MeV: a 3× gap, robust.
rng = np.random.default_rng(0)
al_in = al_tr >= 108; x17_in = x17_tr >= 108
RA_win = al_in.sum() * wA
RX_win = x17_in.sum() * wX
# total pair KE (truth): Al from reco entry KE; signal ~ full transition (fixed)
al_tot = (a["em_ke_entry"] + a["ep_ke_entry"])[al_pair][al_in]
al_tot = al_tot[np.isfinite(al_tot)]
sog = np.load("analysis/al_pair/signal_openingangle.npz")
x17_tot = (sog["x17_em"] + sog["x17_ep"])   # ~19.56 MeV, kinematically fixed
al_emin_win = al_emin[al_in]; x17_emin_win = x17_emin[x17_in]

def frac_above(e, thr, sig, n=60):
    if len(e) == 0: return 0.0
    if sig == 0: return float((e > thr).mean())
    return float(np.mean([(e*(1+sig*rng.standard_normal(len(e))) > thr).mean()
                          for _ in range(n)]))

print("\n== residual Al under X17 window (θ>108°) vs MM momentum resolution ==")
print(f"{'sigma_p/p':>9} | {'softer>4MeV':>22} | {'TOTAL>13MeV':>22}")
print(f"{'':>9} | {'Al/day':>10} {'X17eff':>7} {'Al:X17':>3} | {'Al/day':>10} {'X17eff':>7} {'Al:X17':>3}")
scan = {}
for sig in (0.0, 0.10, 0.20, 0.30, 0.50):
    al_s = frac_above(al_emin_win, 4.0, sig) * RA_win
    x_s  = frac_above(x17_emin_win, 4.0, sig)
    al_t = frac_above(al_tot, 13.0, sig) * RA_win
    x_t  = frac_above(x17_tot, 13.0, sig)
    scan[sig] = (al_s, x_s, al_t, x_t)
    r_s = al_s/(x_s*RX_win) if x_s else np.inf
    r_t = al_t/(x_t*RX_win) if x_t else np.inf
    print(f"{sig*100:>8.0f}% | {al_s:>10.3g} {x_s:>7.2f} {r_s:>8.2g} | {al_t:>10.3g} {x_t:>7.2f} {r_t:>8.2g}")
print(f"(X17 signal in window = {RX_win:.3g}/day; IPC in window = {R_IPC*ipc_acc*(ipc_tr>=108).mean():.3g}/day)")

# analytic Gaussian-tail leak for the TOTAL-energy cut (reaches below the MC
# sample floor ~{RA_win/max(al_in.sum(),1):.0f}/day). Per pair: total_meas ~
# N(E1+E2, sig*sqrt(E1^2+E2^2)); P(>13) via erfc. Sum × per-pair weight.
try:
    from scipy.special import erfc
    e1 = a["em_ke_entry"][al_pair][al_in]; e2 = a["ep_ke_entry"][al_pair][al_in]
    ok = np.isfinite(e1) & np.isfinite(e2)
    e1, e2 = e1[ok], e2[ok]
    tot = e1 + e2; sd1 = np.sqrt(e1**2 + e2**2)
    fine = np.arange(0.10, 0.51, 0.02)
    print("\n== analytic TOTAL>13 MeV leak (Gaussian per-lepton resolution) ==")
    al_ana = []
    for sig in fine:
        z = (13.0 - tot) / (sig * sd1)
        p = 0.5 * erfc(z / np.sqrt(2))
        al_day = float(p.sum()) * wA
        al_ana.append(al_day)
        if abs(sig*100 - round(sig*100/5)*5) < 1e-6:
            print(f"  σ={sig*100:4.0f}%  Al total-cut residual = {al_day:.3g}/day  "
                  f"({'BELOW' if al_day < RX_win else 'above'} X17 {RX_win:.2g}/day)")
    # σ where Al residual crosses the X17 signal rate
    al_ana = np.array(al_ana)
    below = np.where(al_ana < RX_win)[0]
    sig_ok = fine[below[-1]]*100 if len(below) else None
    if sig_ok:
        print(f"  => Al total-cut residual < X17 signal for σ(p)/p ≲ {sig_ok:.0f}%")
except ImportError:
    fine, al_ana, sig_ok = None, None, None

# ── figure ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 3, figsize=(18, 4.8))
b = np.linspace(0, 180, 91); bw = b[1]-b[0]
for tr, w, col, lab in [(al_tr, wA, "0.25", f"Al (MM pairs, {RA:.2g}/day)"),
                        (ipc_tr, wI, "C0", f"IPC ({RI:.2g}/day)"),
                        (x17_tr, wX, "C3", f"X17 ({RX:.2g}/day)")]:
    h, _ = np.histogram(tr, bins=b)
    ax[0].step(b[:-1], h*w, where="post", color=col, lw=2, label=lab)
ax[0].axvline(108, color="C3", ls=":", lw=1.2)
ax[0].set_yscale("log"); ax[0].set_ylim(1e-4, 1e6)
ax[0].set(xlabel="RECONSTRUCTED opening angle [deg]", ylabel="MM pairs / day / bin",
          title="(A) Reconstructed (MSC) opening angle, rate-weighted")
ax[0].legend(fontsize=8.5)

# Panel B: soft-lepton handle on the MM pairs in the X17 window
inwin_al = (al_tr >= 108)
inwin_x = (x17_tr >= 108)
kb = np.linspace(0, 12, 61)
ax[1].hist(al_emin[inwin_al], bins=kb, weights=np.full(inwin_al.sum(), wA),
           histtype="step", lw=2, color="0.25", label="Al (reco θ>108°)")
ax[1].hist(x17_emin[inwin_x], bins=kb, weights=np.full(inwin_x.sum(), wX),
           histtype="step", lw=2, color="C3", label="X17 (reco θ>108°)")
ax[1].axvline(4.0, color="k", ls="--", lw=1)
ax[1].set_yscale("log")
ax[1].set(xlabel="softer-lepton KE (truth) [MeV]", ylabel="pairs / day / bin",
          title="(B) Soft-lepton handle inside the X17 window")
ax[1].legend(fontsize=9)

# Panel C: residual Al vs momentum resolution, both energy handles
sigs = sorted(scan)
xs = [s*100 for s in sigs]
ax[2].plot(xs, [scan[s][0] for s in sigs], "s-", color="0.55", lw=2,
           label="Al resid (softer>4 MeV)")
if fine is not None:
    ax[2].plot(fine*100, np.clip(al_ana, 1e-4, None), "-", color="0.15", lw=2,
               label="Al resid (total>13 MeV, analytic)")
    if sig_ok:
        ax[2].axvline(sig_ok, color="green", ls=":", lw=1.2)
ax[2].plot(xs, [scan[s][2] for s in sigs], "o", color="0.15", ms=6)
ax[2].axhline(RX_win, color="C3", ls="--", lw=1.5, label=f"X17 signal ({RX_win:.2g}/day)")
ax[2].set_yscale("log"); ax[2].set_ylim(1e-4, 1e5)
ax[2].set(xlabel="MM per-lepton momentum resolution σ(p)/p [%]",
          ylabel="rate in X17 window [/day]",
          title="(C) Residual Al vs momentum resolution")
ax[2].legend(fontsize=8.5)
fig.suptitle("Al pair background — detector-level verdict: opening angle is hopeless; "
             "a tracking total-energy cut is the discriminant",
             y=1.02, fontsize=12)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"analysis/al_pair/convpair_verdict.{ext}", bbox_inches="tight", dpi=140)
print("\nwrote analysis/al_pair/convpair_verdict.png/.pdf")
