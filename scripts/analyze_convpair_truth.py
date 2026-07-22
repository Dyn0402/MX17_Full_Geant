#!/usr/bin/env python3
"""analyze_convpair_truth.py — Task 1b/3/4a: truth-level γ→e+e- conversion
background from the ConvPairTree campaign (neutrons_convpair_2cm).

For every conversion it computes the truth opening angle, softer-lepton KE,
total pair energy, and transverse vertex radius, grouped by conversion volume
and by capture-γ line (gamma_E). It rate-normalizes to pairs/pulse and
pairs/day using the beam flux integral over the sampling window, and ranks the
pair-producing materials by rate × proximity to the signal opening-angle window.

Run on lxplus (EOS + uproot). Writes a small npz to /tmp/dneff for local plots.

    python3 scripts/analyze_convpair_truth.py [MIN_FILE_GB]
"""
import sys, glob, os, json, collections
import numpy as np, uproot

DIR   = "/eos/experiment/ntof/data/x17/full_sim/neutrons_convpair_2cm"
NEV_PER_FILE = 10_000_000
FLUX_WINDOW  = 4_292_400.81     # n/pulse in [0.001,2] eV (flux_n_pulse_NOisolet_100bpd)
PULSES_DAY   = 1.929e4
MINGB = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
OUT   = "/tmp/dneff/convpair_truth.npz"

fs = sorted(glob.glob(f"{DIR}/*.root"))
done = [f for f in fs if os.path.getsize(f) > MINGB * 1e9]
print(f"{len(fs)} files present, {len(done)} complete (>{MINGB} GB)")

br = ["conv_vol", "gamma_E", "openingAngle", "em_ke", "ep_ke",
      "vx", "vy", "vz", "neutron_E_eV", "weight", "capture_vol"]
acc = {k: [] for k in br if k not in ("conv_vol", "capture_vol")}
convvol = []; capvol = []
n_events = 0
for fp in done:
    try:
        f = uproot.open(fp)
        a = f["ConvPairTree"].arrays(br, library="np")
        n_events += f["EventTree"].num_entries
    except Exception as e:
        print("  skip", fp.split("/")[-1], repr(e)[:60]); continue
    for k in acc:
        acc[k].append(a[k])
    convvol.append(np.array([x.decode() if isinstance(x, bytes) else x for x in a["conv_vol"]]))
    capvol.append(np.array([x.decode() if isinstance(x, bytes) else x for x in a["capture_vol"]]))
if n_events == 0:
    print("no complete files yet — rerun when the campaign has output"); sys.exit(0)

for k in acc:
    acc[k] = np.concatenate(acc[k])
convvol = np.concatenate(convvol); capvol = np.concatenate(capvol)
w_pulse = FLUX_WINDOW / n_events          # per-sim-neutron → per-pulse weight
theta = acc["openingAngle"]
emin  = np.minimum(acc["em_ke"], acc["ep_ke"])   # softer lepton
etot  = acc["em_ke"] + acc["ep_ke"] + 1.022
rvert = np.hypot(acc["vx"], acc["vz"])           # transverse radius (beam=+Y) [mm]

print(f"\nsim events={n_events:,}  total conv pairs={len(theta):,}  "
      f"weight/pulse={w_pulse:.3e}")
print(f"ALL conversions rate: {len(theta)*w_pulse:.3g}/pulse = "
      f"{len(theta)*w_pulse*PULSES_DAY:.3g}/day")

# ── group by conversion volume ──────────────────────────────────────────────
print("\n{:<18} {:>8} {:>11} {:>10} {:>10} {:>9} {:>9}".format(
    "conv_vol", "N", "/day", "med θ", ">100° /day", "medEγ", "med θ<40"))
rows = []
for vol, n in collections.Counter(convvol).most_common(15):
    m = convvol == vol
    per_day = n * w_pulse * PULSES_DAY
    tail_day = int((theta[m] > 100).sum()) * w_pulse * PULSES_DAY
    rows.append((vol, n, per_day, float(np.median(theta[m])),
                 tail_day, float(np.median(acc["gamma_E"][m]))))
    print("{:<18} {:>8d} {:>11.3g} {:>9.1f}° {:>10.3g} {:>8.2f} {:>8.1f}%".format(
        vol, n, per_day, np.median(theta[m]), tail_day,
        np.median(acc["gamma_E"][m]), (theta[m] < 40).mean() * 100))

# ── Al capsule detail (the dominant background) ─────────────────────────────
al = convvol == "He3Cap_Al"
print(f"\n=== He3Cap_Al ({al.sum()} pairs, {al.sum()*w_pulse*PULSES_DAY:.3g}/day) ===")
print("  gamma_E lines:", collections.Counter(np.round(acc['gamma_E'][al], 2)).most_common(6))
print(f"  truth opening: median={np.median(theta[al]):.1f}  "
      f">100°={(theta[al]>100).mean()*100:.3f}%  >108°(X17)={(theta[al]>108).mean()*100:.3f}%")
print(f"  softer lepton KE: median={np.median(emin[al]):.2f} MeV  "
      f">4MeV(X17 floor)={(emin[al]>4).mean()*100:.3f}%")
print(f"  vertex radius: median={np.median(rvert[al]):.1f} mm  "
      f"[{np.percentile(rvert[al],5):.1f},{np.percentile(rvert[al],95):.1f}]")

np.savez(OUT, theta=theta, emin=emin, etot=etot, rvert=rvert,
         gamma_E=acc["gamma_E"], convvol=convvol, capvol=capvol,
         em_ke=acc["em_ke"], ep_ke=acc["ep_ke"],
         vx=acc["vx"], vy=acc["vy"], vz=acc["vz"],
         n_events=n_events, w_pulse=w_pulse, pulses_day=PULSES_DAY)
meta = {"n_events": n_events, "n_files": len(done), "w_pulse": w_pulse,
        "pulses_day": PULSES_DAY, "flux_window_npulse": FLUX_WINDOW,
        "al_pairs_per_day": float(al.sum() * w_pulse * PULSES_DAY),
        "al_pairs_per_pulse": float(al.sum() * w_pulse)}
with open("/tmp/dneff/convpair_truth.json", "w") as f:
    json.dump(meta, f, indent=2)
print(f"\nwrote {OUT} + convpair_truth.json")
