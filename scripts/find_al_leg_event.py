#!/usr/bin/env python3
"""find_al_leg_event.py — scan a trajdump CSV for an Al-capture event whose
capture gamma deposits in both a SiPM bar (PlasticScint) and a plastic bar
(BackScint) in the SAME arm (a 'leg'). Writes that event's rows to a small CSV.
Arm is inferred from hit position (+X/-X/+Z/-Z)."""
import sys, csv
from collections import defaultdict

FN = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else "/tmp/dneff/al_leg_event.csv"

rows_by_ev = defaultdict(list)
has_alcap = set()
# per event/arm energy in SiPM and plastic
esipm = defaultdict(lambda: defaultdict(float))
eplas = defaultdict(lambda: defaultdict(float))

def arm_of(x, z):
    return ("+X" if x > 0 else "-X") if abs(x) >= abs(z) else ("+Z" if z > 0 else "-Z")

with open(FN) as f:
    r = csv.reader(f); next(r)
    for row in r:
        if len(row) < 16:
            continue
        ev = int(row[0]); vol = row[13]; proc = row[14]
        try:
            edep = float(row[15]); px = float(row[5]); pz = float(row[7])
        except ValueError:
            continue
        rows_by_ev[ev].append(row)
        if vol == "He3Cap_Al" and proc == "nCapture":
            has_alcap.add(ev)
        if edep > 0:
            a = arm_of(px, pz)
            if vol == "PlasticScint":
                esipm[ev][a] += edep
            elif vol in ("BackScintL", "BackScintR"):
                eplas[ev][a] += edep

# score events: Al capture + both SiPM and plastic deposit in same arm
best = None
for ev in has_alcap:
    for a in set(esipm[ev]) & set(eplas[ev]):
        score = min(esipm[ev][a], eplas[ev][a])
        if best is None or score > best[0]:
            best = (score, ev, a, esipm[ev][a], eplas[ev][a])

if best is None:
    print("no Al-capture leg event found")
    sys.exit(1)
score, ev, a, es, ep = best
print(f"best Al-leg event: eventID={ev} arm={a}  "
      f"SiPM edep={es*1000:.1f} keV  plastic edep={ep*1000:.2f} MeV  (min={score:.4f} MeV)")
with open(OUT, "w") as f:
    f.write("eventID,trackID,parentID,particle,istep,pre_x,pre_y,pre_z,"
            "post_x,post_y,post_z,ke_pre_MeV,ke_post_MeV,volume,process,edep_MeV\n")
    for row in rows_by_ev[ev]:
        f.write(",".join(row) + "\n")
print(f"wrote {OUT}  ({len(rows_by_ev[ev])} steps)")
# quick mechanism summary: gamma interaction processes in this event
gproc = defaultdict(int)
for row in rows_by_ev[ev]:
    if row[3] == "gamma" and row[14] not in ("Transportation", "none", ""):
        gproc[row[14]] += 1
print("gamma interaction processes in event:", dict(gproc))
