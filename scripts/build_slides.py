#!/usr/bin/env python3
"""build_slides.py — assemble the MX17 results deck (16:9 .pptx).

Embeds the slide figures (scripts/make_slides_figures.py) plus selected
existing PNGs from docs/. Upload the resulting .pptx to Google Drive with
conversion to Google Slides for the collaborator.

Deck layout
  Section A (INFN request): title, setup, production/reach, July optimistic,
    post-LS3.
  Section B (broader results, for D.N.): production budget, angular-resolution
    floor, the "money plot", resolution vs pair kinematics, the E0 channel,
    backgrounds/timing, summary.

Caveat boxes are added as SEPARATE text shapes (clearly labelled) so the
collaborator can delete them before the INFN talk.

Usage:  /home/dylan/PycharmProjects/nTof_x17/venv/bin/python scripts/build_slides.py
"""

from pathlib import Path

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from PIL import Image

REPO = Path(__file__).resolve().parent.parent
SF = REPO / "docs/slides/figs"           # generated figures
AR = REPO / "docs/angular_resolution/figs"
RP = REPO / "docs/report/figs"
NH = REPO / "docs/neutron_histories/figs"
GEO = REPO / "scripts"
OUT = REPO / "docs/slides/MX17_results.pptx"

# 16:9
SW, SH = Inches(13.333), Inches(7.5)
NAVY = RGBColor(0x1F, 0x3A, 0x5F)
GREY = RGBColor(0x55, 0x55, 0x55)
RED = RGBColor(0xB0, 0x20, 0x20)
CAVEAT_BG = RGBColor(0xFF, 0xF4, 0xD6)
CAVEAT_BORDER = RGBColor(0xD0, 0xA0, 0x30)

prs = Presentation()
prs.slide_width = SW
prs.slide_height = SH
BLANK = prs.slide_layouts[6]


def _txbox(slide, l, t, w, h, lines, *, size=18, bold=False, color=GREY,
           align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, line_spacing=1.05):
    """lines: str or list of (text, dict-overrides) or plain str."""
    box = slide.shapes.add_textbox(l, t, w, h)
    tf = box.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    if isinstance(lines, str):
        lines = [lines]
    for i, ln in enumerate(lines):
        ov = {}
        if isinstance(ln, tuple):
            ln, ov = ln
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = ov.get("align", align)
        p.line_spacing = line_spacing
        if ov.get("space_before"):
            p.space_before = Pt(ov["space_before"])
        r = p.add_run()
        r.text = ln
        f = r.font
        f.size = Pt(ov.get("size", size))
        f.bold = ov.get("bold", bold)
        f.color.rgb = ov.get("color", color)
        if ov.get("italic"):
            f.italic = True
    return box


def _title(slide, text, sub=None):
    bar = slide.shapes.add_shape(1, 0, 0, SW, Inches(0.12))
    bar.fill.solid(); bar.fill.fore_color.rgb = NAVY; bar.line.fill.background()
    _txbox(slide, Inches(0.5), Inches(0.22), Inches(12.9), Inches(0.7),
           text, size=26, bold=True, color=NAVY)
    if sub:
        _txbox(slide, Inches(0.5), Inches(0.95), Inches(12.9), Inches(0.5),
               sub, size=15, color=GREY)


def _img_fit(slide, path, l, t, max_w, max_h):
    """Add image scaled to fit (max_w,max_h) box, centred in it."""
    with Image.open(path) as im:
        iw, ih = im.size
    ar = iw / ih
    box_ar = max_w / max_h
    if ar > box_ar:
        w = max_w; h = Emu(int(max_w / ar))
    else:
        h = max_h; w = Emu(int(max_h * ar))
    l2 = Emu(int(l + (max_w - w) / 2))
    t2 = Emu(int(t + (max_h - h) / 2))
    return slide.shapes.add_picture(str(path), l2, t2, w, h)


def _caveat(slide, l, t, w, h, lines, header="Caveats / assumptions "
            "(separate box — delete before INFN if desired)"):
    box = slide.shapes.add_shape(1, l, t, w, h)
    box.fill.solid(); box.fill.fore_color.rgb = CAVEAT_BG
    box.line.color.rgb = CAVEAT_BORDER; box.line.width = Pt(1)
    tf = box.text_frame; tf.word_wrap = True
    tf.margin_left = Pt(8); tf.margin_right = Pt(8); tf.margin_top = Pt(5)
    p0 = tf.paragraphs[0]; r = p0.add_run(); r.text = header
    r.font.size = Pt(11); r.font.bold = True
    r.font.color.rgb = RGBColor(0x8A, 0x6D, 0x10)
    for ln in lines:
        p = tf.add_paragraph(); p.line_spacing = 1.0
        rr = p.add_run(); rr.text = "•  " + ln
        rr.font.size = Pt(11); rr.font.color.rgb = RGBColor(0x6B, 0x55, 0x10)
    return box


def _notes(slide, text):
    slide.notes_slide.notes_text_frame.text = text


def add(title=None, sub=None):
    s = prs.slides.add_slide(BLANK)
    if title:
        _title(s, title, sub)
    return s


# ════════════════════════════════════════════════════════════════════════════
# SECTION A — the INFN request
# ════════════════════════════════════════════════════════════════════════════

# 1 — TITLE -------------------------------------------------------------------
s = add()
band = s.shapes.add_shape(1, 0, Inches(2.2), SW, Inches(2.9))
band.fill.solid(); band.fill.fore_color.rgb = NAVY; band.line.fill.background()
_txbox(s, Inches(0.8), Inches(2.5), Inches(11.7), Inches(1.4),
       "The X17 search at n_TOF EAR2", size=42, bold=True,
       color=RGBColor(0xFF, 0xFF, 0xFF))
_txbox(s, Inches(0.8), Inches(3.7), Inches(11.7), Inches(1.0),
       "Full Geant4 simulation: production, acceptance, and the path to a "
       "measurement\nFinal July setup  ·  expected results  ·  the post-LS3 "
       "upgrade", size=19, color=RGBColor(0xD7, 0xE3, 0xF4))
_txbox(s, Inches(0.8), Inches(5.4), Inches(11.7), Inches(0.6),
       "D. Neff — 16 June 2026   |   simulation + analysis: this repository "
       "(MX17_Full_Geant)", size=14, color=GREY)
_notes(s, "Deck for INFN presentation of the July run plan + post-LS3 "
       "expectations, plus a broader set of simulation results. Section A "
       "(slides 2-5) is the INFN-requested core; Section B is supporting "
       "detail. All numbers from the full Geant4 campaign (5e8 neutrons, "
       "1e7 pair events).")

# 2 — SETUP -------------------------------------------------------------------
s = add("The experiment: a 4-arm e⁺ e⁻ pair spectrometer",
        "³He(n,γ)⁴He* → (X17 or IPC) → e⁺ e⁻ ; "
        "measure the opening angle in 4 Micromegas + scintillator arms")
_img_fit(s, GEO / "mx17_geometry_topdown.png",
         Inches(0.3), Inches(1.55), Inches(7.4), Inches(5.7))
_txbox(s, Inches(7.9), Inches(1.7), Inches(5.1), Inches(5.3), [
    ("Target", {"size": 18, "bold": True, "color": NAVY}),
    "³He gas capsule, 500 bar, r = 10 mm, 60 mm on axis; thin Al (0.6 mm) "
    "+ CFRP (0.9 mm) wall.",
    ("Tracking", {"size": 18, "bold": True, "color": NAVY,
                  "space_before": 8}),
    "4 Micromegas arms, front face ≈ 22 cm from the target, 30 mm drift "
    "gap — the e⁺ e⁻ directions.",
    ("Calorimetry / trigger", {"size": 18, "bold": True, "color": NAVY,
                               "space_before": 8}),
    "Liquid-scintillator + plastic layers per arm (assumed saturated/unusable "
    "for the high-energy measurement — tracking-only baseline).",
    ("The observable", {"size": 18, "bold": True, "color": NAVY,
                        "space_before": 8}),
    "X17 (m ≈ 17 MeV) makes a sharp opening-angle shoulder at θ ≈ "
    "109°; the IPC continuum is the irreducible background.",
], size=15, color=GREY)
_notes(s, "Geometry from scripts/plot_geometry.py (SimConfig defaults). Beam "
       "into the page. Production happens in the He-3 gas; the two leptons exit "
       "through the capsule wall (the dominant resolution limit, see later) and "
       "are tracked in the Micromegas. Calorimetry is conservatively treated as "
       "unusable at high energy, so the measurement rests on MM tracking only.")

# 3 — PRODUCTION / REACH ------------------------------------------------------
s = add("Production is settled — and it lives at MeV energies",
        "Directly counted ³He(n,γ) from 5×10⁸ neutrons; "
        "the γ-flash sets the high-energy reach")
_img_fit(s, SF / "fig_production_windows.png",
         Inches(0.3), Inches(1.55), Inches(8.2), Inches(5.7))
_txbox(s, Inches(8.7), Inches(1.7), Inches(4.4), Inches(5.4), [
    ("Sub-keV is dead.", {"size": 16, "bold": True, "color": RED}),
    "The gas is opaque to (n,p); radiative capture is capped at "
    "σₙγ/σₙₚ ≈ 10⁻⁸ → ~0.1 "
    "X17/day.",
    ("The MeV window is the measurement.", {"size": 16, "bold": True,
                                            "color": NAVY, "space_before": 10}),
    "Above ~100 keV the gas is thin and the rate table is confirmed by Geant4 "
    "(within 10–70%). ~33 X17/day produced over the full range.",
    ("The γ-flash limits the reach.", {"size": 16, "bold": True,
                                            "color": NAVY, "space_before": 10}),
    "Faster neutrons arrive sooner, closer to the flash. How early the "
    "detectors recover sets the max neutron energy — this is the LS3 "
    "lever.",
], size=14, color=GREY)
_notes(s, "From docs/report/mev_note.tex + analysis/mev. The per-bin X17/day "
       "uses alpha_IPC = 3.5e-3. Green = July reach (0.2-0.7 MeV), blue = the "
       "extra reach LS3 buys (0.7-2 MeV). Top axis is neutron TOF to EAR2 at "
       "19.5 m: 2 MeV = 1.0 us, 0.7 MeV = 1.7 us, 0.2 MeV = 3.2 us. Recorded "
       "X17 in 30 d: 64 (July) -> 220 (LS3).")

# 4 — JULY OPTIMISTIC ---------------------------------------------------------
s = add("July run: expected result (optimistic scenario)",
        "Window 0.2–0.7 MeV — conservative γ-flash recovery "
        "(detectors usable by ≈ 1.7 µs); 30-day run")
_img_fit(s, SF / "fig_stacked_july.png",
         Inches(0.25), Inches(1.5), Inches(8.0), Inches(4.6))
_txbox(s, Inches(8.4), Inches(1.6), Inches(4.7), Inches(3.3), [
    ("≈ 64 X17 on ≈ 3 100 IPC", {"size": 20, "bold": True,
                                           "color": NAVY}),
    "recorded in 30 days (S/B ≈ 0.021).",
    ("The peak is diluted, not visible by eye.", {"size": 15, "bold": True,
                                                  "color": RED,
                                                  "space_before": 8}),
    "Capsule multiple scattering smears θ by σ68 ≈ 14°; the "
    "109° shoulder (green dotted = truth) becomes a broad bump. The "
    "measurement is a template fit, not a cut-and-count.",
], size=14, color=GREY)
_caveat(s, Inches(0.4), Inches(6.25), Inches(12.5), Inches(1.05), [
    "Optimistic: best-estimator (target-centre chord) smearing; no pile-up, "
    "no γ-flash losses, no θ-dependent MM acceptance shaping.",
    "α_IPC = 3.5×10⁻³ (¹²C/⁸Be-anchored) and "
    "X17/IPC = 2.5% are assumed inputs (provenance pending). Table α = "
    "2.1×10⁻³ would scale signal ×0.6 (→ ~38 X17).",
    "At-rest kinematics (shoulder shifts a few degrees with the Eₙ boost). "
    "MM-double acceptance 19.6% (X17) / 23.6% (IPC).",
])
_notes(s, "Conservative July choice: gamma-flash recovery only by ~1.7 us, so "
       "max neutron energy ~0.7 MeV (this 'energy cutoff' framing is how the "
       "collaboration thinks). 143 captures in the window -> 64 recorded X17 / "
       "3093 IPC over 30 days at alpha=3.5e-3. The caveat box is a separate "
       "shape: delete it for the INFN talk if desired. Figure: "
       "scripts/make_slides_figures.py.")

# 5 — POST-LS3 ----------------------------------------------------------------
s = add("After LS3: extend the reach to 2 MeV neutrons",
        "Faster γ-flash recovery (≈ 1.0 µs) → window 0.2–2 MeV; "
        "same analysis, ×3.4 the statistics")
_img_fit(s, SF / "fig_stacked_compare.png",
         Inches(0.2), Inches(1.6), Inches(12.9), Inches(4.55))
_txbox(s, Inches(0.5), Inches(6.2), Inches(12.4), Inches(1.1), [
    ("July 64 X17 → post-LS3 220 X17 recorded (30 d), on the same "
     "S/B ≈ 0.021.", {"size": 17, "bold": True, "color": NAVY}),
    "Reaching 2 MeV neutrons (faster γ-flash recovery) adds the most "
    "productive energy decade. The spectrum shape is unchanged — the gain "
    "is purely statistical, which is exactly what a template fit converts into "
    "significance.",
], size=14, color=GREY)
_notes(s, "LS3 goal stated by the collaboration: achieve neutron energy of at "
       "least 2 MeV. Same plot as July, window extended 0.7->2 MeV. Recorded "
       "X17 30 d: 220 vs 64 (x3.4). S/B identical (set by X17/IPC branching x "
       "acceptance ratio, window-independent). The decisive question remains "
       "the gamma-flash recovery time, measured on the real setup.")

# ════════════════════════════════════════════════════════════════════════════
# SECTION B — broader results
# ════════════════════════════════════════════════════════════════════════════

# 6 — SECTION DIVIDER ---------------------------------------------------------
s = add()
band = s.shapes.add_shape(1, 0, Inches(3.0), SW, Inches(1.5))
band.fill.solid(); band.fill.fore_color.rgb = NAVY; band.line.fill.background()
_txbox(s, Inches(0.8), Inches(3.25), Inches(11.7), Inches(1.0),
       "Supporting results", size=36, bold=True,
       color=RGBColor(0xFF, 0xFF, 0xFF))
_txbox(s, Inches(0.8), Inches(4.7), Inches(11.7), Inches(0.6),
       "Production budget · the angular-resolution floor · the analysis "
       "strategy · the γ-dark E0 channel · backgrounds", size=16,
       color=GREY)

# 6b — GEOMETRY: SIDE VIEW + 3D ----------------------------------------------
s = add("Detector geometry — side view and 3D",
        "Companion views to the top-down layout: the ³He capsule shape, "
        "beam-end-on incidence, and the full 4-arm arrangement")
_img_fit(s, GEO / "mx17_geometry_sideview.png",
         Inches(0.25), Inches(1.55), Inches(7.1), Inches(5.4))
_img_fit(s, GEO / "mx17_geometry_3d.png",
         Inches(7.5), Inches(1.55), Inches(5.6), Inches(5.4))
_txbox(s, Inches(0.4), Inches(7.0), Inches(6.7), Inches(0.4),
       "Side view (XY): beam enters end-on along +Y into the capsule; "
       "arms 0/1 (±X) shown.", size=12, color=GREY)
_txbox(s, Inches(7.5), Inches(7.0), Inches(5.6), Inches(0.4),
       "3D: 4 arms (±X, ±Z), each MM + scintillator stack around the target.",
       size=12, color=GREY)
_notes(s, "Backup geometry slide (scripts/plot_geometry.py). Side view shows "
       "the capsule bottle profile (20 mm dia x 40 mm + domed ends, neck/valve "
       "on axis) and the beam hitting it end-on along +Y. The 3D isometric "
       "shows all four arms; arms 2,3 are along +/-Z (not shown in the side "
       "view).")

# 7 — ANGULAR RESOLUTION FLOOR ------------------------------------------------
s = add("The angular resolution is hardware-frozen",
        "Opening-angle resolution is set by multiple scattering in the "
        "capsule wall — no reconstruction can beat it")
_img_fit(s, AR / "fig_ms_budget.png",
         Inches(0.3), Inches(1.6), Inches(7.6), Inches(5.5))
_txbox(s, Inches(8.1), Inches(1.7), Inches(4.9), Inches(5.4), [
    ("σ68(Δθ) ≈ 14.5°", {"size": 26, "bold": True,
                                                 "color": NAVY}),
    "full accepted spectrum (≈ 11.5–12.5° for symmetric pairs), "
    "even assuming a perfect vertex.",
    ("Why:", {"size": 16, "bold": True, "color": NAVY, "space_before": 10}),
    "The 0.6 mm Al + 0.9 mm CFRP wall (~11 mm from the vertex) is 79% of the "
    "material budget. A scatter that early rotates the whole visible track.",
    ("Verdict:", {"size": 16, "bold": True, "color": RED, "space_before": 10}),
    "All 7 reconstruction methods land within 2° of each other. The only "
    "hardware lever is an all-composite vessel (→ ~10°). The recovery "
    "path is statistical: fit the shape.",
], size=14, color=GREY)
_notes(s, "From docs/angular_resolution/angular_resolution_note. Highland "
       "budget on the new STEP capsule. The vertex-constraint trick "
       "(target-centre chord) is already within 0.5 deg of the perfect-vertex "
       "oracle because the small target banks the gain. Detector spatial "
       "resolution costs <0.5 deg. This is the single most important hardware "
       "fact for the measurement.")

# 8 — THE MONEY PLOT ----------------------------------------------------------
s = add("Why the analysis is a template fit, not a peak cut",
        "Truth → stacked on IPC → both smeared: the shoulder washes "
        "out, but its shape is known to ~1% from 10⁷ simulated pairs")
_img_fit(s, AR / "fig_theta_money.png",
         Inches(0.25), Inches(1.6), Inches(12.85), Inches(4.7))
_txbox(s, Inches(0.5), Inches(6.35), Inches(12.4), Inches(0.95), [
    "Cut-and-count on θ is dead (recorded S/B ≈ 0.02). What survives: "
    "a fit of the known smeared-signal template + the smooth IPC continuum, "
    "with resolution entering as dilution of statistical power rather than as "
    "a cliff. This is the natural go/no-go calculation.",
], size=14, color=GREY)
_notes(s, "The 'money plot' from the angular note. Per-leg smearing drawn from "
       "the measured P(psi|KE) response (best estimator). Normalised to "
       "full-production 30-day recorded yields (191 X17 / 9190 IPC at "
       "alpha=2.1e-3). The middle vs right panels show the shoulder before/after "
       "smearing on the IPC continuum.")

# 9 — RESOLUTION VS KINEMATICS ------------------------------------------------
s = add("Resolution is a soft-track problem",
        "σ68 vs pair symmetry, and the +bias (correctable in a "
        "response-matrix fit) — X17 is sharpest at its symmetric 109° shoulder")
_img_fit(s, AR / "fig_sigma68_vs_minke.png",
         Inches(0.4), Inches(1.7), Inches(12.5), Inches(5.0))
_txbox(s, Inches(0.5), Inches(6.75), Inches(12.4), Inches(0.55),
       "Below ~3 MeV a leg is directionless after the wall (σ68 → "
       "40°); above ~8 MeV it plateaus at ~12°. X17 is sharpest at "
       "its symmetric 109° shoulder — the good news.",
       size=14, color=GREY)
_notes(s, "X17 has KE- + KE+ ~ 19.6 MeV so one leg is always <= 9.8 MeV; the "
       "resolution is driven by the softer leg. MM-internal quality cuts "
       "(straightness, dE/dx) could reject soft legs and recover ~12.5 deg "
       "without a calorimeter -- the one open reconstruction study.")

# 10 — E0 CHANNEL -------------------------------------------------------------
s = add("The γ-dark E0 channel: completeness, not a game-changer",
        "0⁺→0⁺ captures make no photon but can make a massive "
        "X17 — yet the rate is tiny and sub-keV")
_img_fit(s, RP / "fig_e0_final_metric.png",
         Inches(0.3), Inches(1.6), Inches(8.3), Inches(5.4))
_txbox(s, Inches(8.7), Inches(1.7), Inches(4.4), Inches(5.4), [
    ("E0 is real but small.", {"size": 16, "bold": True, "color": NAVY}),
    "The 0⁺ doorway sits sub-keV, where the same 10⁻⁸ radiative "
    "branching that kills the thermal program applies — times f ≤ "
    "10⁻².",
    ("The flash readout cannot see it.", {"size": 16, "bold": True,
                                          "color": RED, "space_before": 10}),
    "Sub-keV → TOF > 45 µs → arrives after the window closes. "
    "Only a dedicated thermal trigger records it (≤ 0.06 X17/day).",
    ("So:", {"size": 16, "bold": True, "color": NAVY, "space_before": 10}),
    "Leaving E0 out does not bias the baseline high-energy analysis. Worth "
    "adding to the generator for completeness only.",
], size=14, color=GREY)
_notes(s, "From docs/report/e0_final_metric_note + e0_branch docs. Addresses "
       "Alberto's E0 question. Key physics: 0+->0+ is photon-dark because the "
       "photon is massless, but a massive vector/scalar X17 is allowed. The "
       "rate is bounded by the (tiny) E0 capture rate. Left panel: IPC "
       "background pairs/day; right: X17 signal. Gold/blue = thermal vs flash "
       "trigger regions.")

# 11 — SUMMARY ----------------------------------------------------------------
s = add("Summary & status")
_txbox(s, Inches(0.7), Inches(1.6), Inches(12.0), Inches(5.4), [
    ("Production is settled.", {"size": 20, "bold": True, "color": NAVY}),
    "Full Geant4 (5×10⁸ neutrons) confirms the rate table where it is "
    "valid: ~33 X17/day produced, living at MeV energies. The sub-keV "
    "measurement is dead (self-shielding).",
    ("July (current hardware): ~64 recorded X17 / 30 d.", {"size": 20,
        "bold": True, "color": NAVY, "space_before": 12}),
    "Window 0.2–0.7 MeV under a conservative γ-flash recovery; on "
    "~3 100 IPC (S/B ≈ 0.02). Extractable by a template fit, not a peak "
    "cut.",
    ("Post-LS3: ~220 recorded X17 / 30 d.", {"size": 20, "bold": True,
        "color": NAVY, "space_before": 12}),
    "Reaching 2 MeV neutrons (×3.4 statistics). The decisive measurement "
    "is the real γ-flash recovery time.",
    ("The resolution is hardware-frozen at σ68 ≈ 14°.",
     {"size": 20, "bold": True, "color": NAVY, "space_before": 12}),
    "Capsule wall MS; an all-composite vessel is the only lever. E0 is a "
    "completeness item, not a driver.",
    ("Open items:", {"size": 16, "bold": True, "color": RED,
                     "space_before": 12}),
    "α_IPC & X17/IPC provenance (Alberto); Eₙ-dependent generator "
    "kinematics; the template-fit significance projection; γ-flash "
    "recovery measurement.",
], size=15, color=GREY)
_notes(s, "Wrap-up. The program is now an acceptance-and-backgrounds problem, "
       "not a production one. Next analysis step: wire the smeared templates "
       "into the fast-MC likelihood fit and quote expected significance vs run "
       "time.")

prs.save(str(OUT))
print(f"wrote {OUT}  ({len(prs.slides.__iter__.__self__._sldIdLst)} slides)")
