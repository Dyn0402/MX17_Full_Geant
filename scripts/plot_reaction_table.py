#!/usr/bin/env python3
"""
plot_reaction_table.py — reference table of the neutron reactions in the MX17
detector and the particles/energies they produce.  Thermal-capture Q-values are
the neutron separation energy of the product (AME2020 / ENDF-B/VIII.0); the
capture gamma is emitted as a cascade summing to Q.  Signal and the dominant
background are highlighted.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent.parent / "analysis/thermal_2cm"

# columns: Material | Reaction | Target (abund) | Q [MeV] | Products (energy) |
#          gamma emitted | Role
ROWS = [
    ("He-3 gas", "$^{3}$He(n,p)$^{3}$H", "$^{3}$He", "0.764",
     "p 573 keV + t 191 keV", "none (charged)", "Beam monitor (dominant gas rxn)"),
    ("He-3 gas", "$^{3}$He(n,$\\gamma$)$^{4}$He", "$^{3}$He", "20.578",
     "20.58 MeV $\\gamma$  or  e$^+$e$^-$ pair", "20.58 MeV (or IPC)", "SIGNAL — X17 / IPC channel"),
    ("Al capsule", "$^{27}$Al(n,$\\gamma$)$^{28}$Al", "$^{27}$Al (100%)", "7.725",
     "$\\gamma$ cascade $\\to$ Compton e$^-$", "$\\Sigma$ = 7.72 MeV (1779 keV +…)",
     "MAIN coincidence bkg (96% of legs)"),
    ("LS / plastics / CFRP", "$^{1}$H(n,$\\gamma$)$^{2}$H", "$^{1}$H", "2.224",
     "single 2.223 MeV $\\gamma$ $\\to$ e$^-$", "2.223 MeV",
     "Soft $\\gamma$: many singles, few legs"),
    ("plastic / LS (fast n)", "$^{1}$H(n,n')p", "$^{1}$H", "elastic",
     "recoil proton (direct)", "none", "Non-$\\gamma$ ionization bkg"),
    ("CFRP / plastic", "$^{12}$C(n,$\\gamma$)$^{13}$C", "$^{12}$C (98.9%)", "4.946",
     "$\\gamma$ cascade $\\to$ e$^-$", "4.945, 3.089, 1.262 MeV", "Minor $\\gamma$ bkg"),
    ("PMT glass", "$^{28}$Si(n,$\\gamma$)$^{29}$Si", "$^{28}$Si (92%)", "8.473",
     "$\\gamma$ cascade $\\to$ e$^-$", "$\\Sigma$ = 8.47 MeV", "Minor bkg"),
    ("readout Cu", "$^{63}$Cu(n,$\\gamma$)$^{64}$Cu", "$^{63}$Cu (69%)", "7.916",
     "$\\gamma$ cascade $\\to$ e$^-$", "$\\Sigma$ = 7.92 MeV", "Minor bkg"),
    ("air (trace)", "$^{14}$N(n,p)$^{14}$C", "$^{14}$N (99.6%)", "0.626",
     "proton 584 keV", "none", "Trace"),
]
COLS = ["Material / volume", "Reaction", "Target (abund.)", "Q [MeV]",
        "Products (energy)", "$\\gamma$ emitted", "Role"]

fig, ax = plt.subplots(figsize=(15.5, 4.6))
ax.axis("off")
tbl = ax.table(cellText=ROWS, colLabels=COLS, cellLoc="left", loc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(9)
tbl.scale(1, 1.9)
# column widths
widths = [0.13, 0.15, 0.11, 0.06, 0.20, 0.16, 0.19]
for (r, c), cell in tbl.get_celld().items():
    cell.set_width(widths[c])
    cell.set_edgecolor("0.8")
    if r == 0:
        cell.set_facecolor("#34495e"); cell.set_text_props(color="white", weight="bold")
    else:
        # highlight signal (row idx 2 in table = data row 1) and main bkg (row 3)
        role = ROWS[r - 1][-1]
        if role.startswith("SIGNAL"):
            cell.set_facecolor("#d6eaf8")
        elif role.startswith("MAIN"):
            cell.set_facecolor("#fdebd0")
        elif r % 2 == 0:
            cell.set_facecolor("#f6f6f6")

ax.set_title("Neutron reactions in the MX17 detector — products & $\\gamma$ energies "
             "(thermal capture Q = $S_n$ of product; AME2020/ENDF)",
             fontsize=12, pad=14)
fig.text(0.5, 0.02, "Detector hits are deposited by the secondary electrons of "
         "the capture $\\gamma$ (Compton/photoelectric), except the direct recoil "
         "proton from fast-neutron elastic scattering and the charged $^3$He(n,p)t "
         "products (confined to the gas).", ha="center", fontsize=8.5, color="0.35")
fig.savefig(OUT / "reaction_table.pdf", bbox_inches="tight")
fig.savefig(OUT / "reaction_table.png", dpi=150, bbox_inches="tight")
print(f"wrote {OUT}/reaction_table.pdf/.png")
