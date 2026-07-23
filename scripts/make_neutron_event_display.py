#!/usr/bin/env python3
"""
make_neutron_event_display.py

Event displays of neutron histories in the He-3 target, from the --trajdump
CSV(s) produced by mx17_full_sim.  Each neutron path is drawn in the (y, x)
plane (beam axis = +y) over the real capsule cross-section, coloured by
log10(kinetic energy) so the slowing-down is visible; charged secondaries
(p, t from (n,p)t; e+/e- from a pair) and the terminal vertex are marked.

Usage:
    python3 scripts/make_neutron_event_display.py <traj.csv> [more.csv ...] \
        -o docs/neutron_histories/figs/event_display.pdf [--max-events 6]

Each input CSV is one job (typically one energy window); the panel title shows
the incident neutron energy read from the primary's first step.  Columns
(written by SteppingAction::DumpTrajectoryStep):
    eventID,trackID,parentID,particle,istep,
    pre_x,pre_y,pre_z,post_x,post_y,post_z,
    ke_pre_MeV,ke_post_MeV,volume,process,edep_MeV   (positions in mm)
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

# ── Capsule cross-section, copied from DetectorConstruction.cc ────────────────
# Polycones are built in a local frame then rotated rotateX(+90°) (nose-first
# mounting, 2026-07-23 audit), so the local profile axis z -> world y (tip at
# y=-35 mm faces the beam) and the profile radius -> world transverse (x,z).
# We draw the y-x section: world_y = z_profile[mm], world_x = ±r_profile[mm].
Z_GAS = np.array([-30, -28, -26, -24, -22, -20, 20, 22, 24, 26, 28, 30.0])
RO_GAS = np.array([0.001, 6.0, 8.0, 9.165, 9.798, 10.0,
                   10.0, 9.798, 9.165, 8.0, 6.0, 0.001])
Z_VESSEL = np.array([-35.000, -34.635, -33.947, -32.880, -31.447, -29.693,
                     -27.688, -25.521, -23.288, -21.075, -20.000, 20.000,
                     22.000, 24.000, 26.000, 26.770, 31.370, 40.900, 50.980])
RO_AL = np.array([0.000, 2.323, 3.902, 5.433, 6.850, 8.090, 9.102, 9.856,
                  10.342, 10.573, 10.600, 10.600, 10.317, 9.469, 8.054, 7.360,
                  6.912, 3.500, 3.500])
RO_CFRP = np.array([0.000, 3.223, 4.802, 6.333, 7.750, 8.990, 10.002, 10.756,
                    11.242, 11.473, 11.500, 11.500, 11.217, 10.369, 8.954,
                    8.260, 7.812, 4.400, 4.400])

SECONDARY_STYLE = {        # particle -> (color, label)
    "proton":  ("tab:red",    "p"),
    "triton":  ("tab:orange", "t"),
    "e-":      ("tab:green",  "e$^-$"),
    "e+":      ("tab:purple", "e$^+$"),
    "alpha":   ("tab:brown",  r"$\alpha$"),
    "gamma":   ("0.5",        r"$\gamma$"),
    "deuteron":("tab:pink",   "d"),
}


def draw_capsule(ax):
    """Outline gas / Al / CFRP in the (y, x) plane."""
    def band(z, ro, **kw):
        y = np.r_[z, z[::-1]]
        x = np.r_[ro, -ro[::-1]]
        ax.fill(y, x, **kw)
    band(Z_VESSEL, RO_CFRP, facecolor="0.15", edgecolor="0.1", lw=0.8,
         alpha=0.55, zorder=0, label="CFRP wrap")
    band(Z_VESSEL, RO_AL, facecolor="0.75", edgecolor="0.4", lw=0.8,
         alpha=0.7, zorder=1, label="Al vessel")
    band(Z_GAS, RO_GAS, facecolor="#bfe3ff", edgecolor="#3a78a0", lw=1.0,
         alpha=0.7, zorder=2, label="$^3$He gas")


def load_events(path):
    """Return {eventID: {trackID: dict(particle, segs[Nx2x2], ke[N], end)}}."""
    events = defaultdict(lambda: defaultdict(
        lambda: {"particle": None, "y0": [], "x0": [], "y1": [], "x1": [],
                 "ke": [], "end_proc": None, "end_yx": None, "cap_yx": None}))
    with open(path) as f:
        for r in csv.DictReader(f):
            eid, tid = int(r["eventID"]), int(r["trackID"])
            t = events[eid][tid]
            t["particle"] = r["particle"]
            t["y0"].append(float(r["pre_y"]));  t["x0"].append(float(r["pre_x"]))
            t["y1"].append(float(r["post_y"])); t["x1"].append(float(r["post_x"]))
            t["ke"].append(float(r["ke_pre_MeV"]))
            proc = r["process"]
            if proc not in ("Transportation", "none", ""):
                t["end_proc"] = proc
                t["end_yx"] = (float(r["post_y"]), float(r["post_x"]))
            # the physically meaningful "excitation" vertex: 4He* formation
            if proc in ("neutronInelastic", "nCapture"):
                t["cap_yx"] = (float(r["post_y"]), float(r["post_x"]))
    return events


def incident_energy_MeV(tracks):
    """Primary neutron (parentID 0) first-step pre-KE."""
    for t in tracks.values():
        if t["particle"] == "neutron" and t["ke"]:
            return max(t["ke"])      # first step carries the launch energy
    return None


def event_interest(tracks):
    """Rank key: richness of the in-gas history (steps whose pre-point is in the
    gas) first — this surfaces multi-scatter slowing-down where it exists (MeV)
    and capture-in-gas where that's all there is (thermal) — then capture as a
    tiebreak."""
    n_ingas = 0
    captured = False
    for t in tracks.values():
        if t["particle"] != "neutron":
            continue
        for x0, y0 in zip(t["x0"], t["y0"]):
            if in_gas(y0, x0):
                n_ingas += 1
        if (t["end_proc"] in ("neutronInelastic", "nCapture") and t["end_yx"]
                and in_gas(*t["end_yx"])):
            captured = True
    return (n_ingas, captured)


def in_gas(y, x):
    """Point inside the gas capsule cross-section (y along axis, x radial)."""
    if abs(y) <= 20:
        return abs(x) <= 10.0
    if abs(y) <= 30:
        return x * x + (abs(y) - 20.0) ** 2 <= 100.0
    return False


# fixed colour scale = fraction of incident KE remaining (shows slowing-down
# the same way in every panel regardless of the incident energy)
KE_NORM = LogNorm(vmin=1e-3, vmax=1.0)


def draw_event(ax, tracks, E_inc):
    draw_capsule(ax)
    E_inc = E_inc or 1.0
    lc_ref = None
    for t in tracks.values():
        segs = np.stack([np.c_[t["y0"], t["x0"]], np.c_[t["y1"], t["x1"]]], axis=1)
        if t["particle"] == "neutron":
            lc = LineCollection(segs, cmap="viridis", norm=KE_NORM, lw=1.8,
                                zorder=5)
            lc.set_array(np.clip(np.array(t["ke"]) / E_inc, 1e-3, 1.0))
            ax.add_collection(lc)
            lc_ref = lc
        else:
            col, _ = SECONDARY_STYLE.get(t["particle"], ("0.3", t["particle"]))
            ax.add_collection(LineCollection(segs, colors=col, lw=1.2,
                                             alpha=0.9, zorder=4))
    # excitation vertices: 4He* formation, i.e. (n,p)t or radiative capture
    for t in tracks.values():
        if t["particle"] == "neutron" and t["cap_yx"]:
            ax.plot(*t["cap_yx"], "*", color="red", ms=14, mec="k", zorder=7)
    ax.set_xlim(-45, 60)
    ax.set_ylim(-22, 22)
    ax.set_xlabel("y  [mm]  (beam →)")
    ax.set_ylabel("x  [mm]")
    ax.set_aspect("equal")
    return lc_ref


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="trajectory CSV file(s)")
    ap.add_argument("-o", "--output", default="event_display.pdf")
    ap.add_argument("--max-events", type=int, default=6,
                    help="events to draw per input file")
    args = ap.parse_args()

    panels = []   # (title, tracks)
    for path in args.inputs:
        events = load_events(path)
        # most interesting events first (gas interaction, then most steps)
        ranked = sorted(events, key=lambda e: event_interest(events[e]),
                        reverse=True)
        for eid in ranked[:args.max_events]:
            tracks = events[eid]
            E = incident_energy_MeV(tracks)
            etxt = (f"{E*1e6:.0f} eV" if E and E < 1e-3 else
                    f"{E*1e3:.1f} keV" if E and E < 1 else
                    f"{E:.2f} MeV" if E else "?")
            nsteps = sum(len(t["ke"]) for t in tracks.values()
                         if t["particle"] == "neutron")
            caps = {t["end_proc"] for t in tracks.values()
                    if t["particle"] == "neutron" and t["cap_yx"]}
            fate = ("$^3$He(n,p)t" if "neutronInelastic" in caps else
                    r"$^3$He(n,$\gamma$)" if "nCapture" in caps else
                    "scattered, escaped")
            panels.append((f"$E_n$={etxt} | {nsteps} n-steps | {fate}",
                           tracks, E))

    if not panels:
        print("No events found in inputs."); return
    n = len(panels)
    ncol = min(2, n)
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2 * ncol, 3.6 * nrow),
                             squeeze=False)
    lc_ref = None
    for ax, (title, tracks, E_inc) in zip(axes.flat, panels):
        lc = draw_event(ax, tracks, E_inc)
        lc_ref = lc_ref or lc
        ax.set_title(title, fontsize=9.5)
    for ax in axes.flat[n:]:
        ax.axis("off")
    if lc_ref is not None:
        cb = fig.colorbar(lc_ref, ax=axes.ravel().tolist(), shrink=0.6,
                          pad=0.02)
        cb.set_label("neutron KE / incident  (slowing-down)")
    # legend proxies
    handles = [plt.Line2D([], [], color=c, lw=2, label=l)
               for c, l in [("#3a78a0", "$^3$He gas"), ("0.75", "Al"),
                            ("0.15", "CFRP")]]
    handles += [plt.Line2D([], [], marker="*", color="red", mec="k", ls="",
                           ms=11, label="$^4$He$^*$ vertex (capture / (n,p)t)")]
    handles += [plt.Line2D([], [], color="tab:red", lw=2, label="p / t / e$^\\pm$")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.02))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", dpi=150)
    if out.suffix != ".png":
        fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved {n} event panel(s) -> {out}")


if __name__ == "__main__":
    main()
