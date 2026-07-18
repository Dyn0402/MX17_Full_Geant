#pragma once
// SimConfig.hh — geometry and run parameters for the 4-arm X17 simulation

#include <string>
#include "G4SystemOfUnits.hh"

struct SimConfig {
    // ── Output ─────────────────────────────────────────────
    std::string outFile  = "x17_output";
    long        seed     = 0;
    int         nEvents  = 10000;
    int         nThreads = 1;
    bool        verbose  = false;

    // ── Generator ──────────────────────────────────────────
    // X17 pair mode
    double x17Mass_MeV          = 16.8;   // X17 boson mass [MeV]
    double transition_energy_MeV = 20.58; // 4He* transition energy [MeV] — used by both X17 and IPC

    // IPC fraction: probability that a given pair event is IPC background (vs. X17 signal).
    // 0 = all X17, 1 = all IPC, 0.5 = equal statistics (recommended for building event pools).
    double ipc_fraction = 0.5;

    // Pair-vertex library: CSV of He3Gas capture positions (volume,x_mm,y_mm,z_mm,
    // from make_capture_library.py --gas-lib).  When set, X17/IPC vertices are
    // sampled from these rows instead of uniformly in the gas — this is how the
    // thermal self-shielding absorption profile (captures within ~mm of the gas
    // entrance face) enters the signal kinematics.  Empty = uniform (default).
    std::string pairVertexLibFile;

    bool        singleParticle            = false;
    std::string singleParticleName        = "e-";
    double      singleParticleEnergy_MeV  = 8.0;
    double      singleParticleTheta_deg   = 90.0;
    double      singleParticlePhi_deg     = 0.0;

    // ── Neutron-beam mode (event_type = 2) ──────────────────
    // Fires neutrons along +Y from y = neutronGunY_cm with energy sampled
    // from the EAR2 evaluated flux and transverse position from the
    // energy-dependent radial profile (both ROOT files in data/).
    bool        neutronMode    = false;
    std::string neutronFluxFile;     // fluxEAR2-Ph3_in_different_units.root
    std::string neutronProfileFile;  // lamda2DvsEn_EAR2.root
    double      neutronEmin_eV = 1e-3;    // sampling window
    double      neutronEmax_eV = 1000.0;  // default: < 1 keV (X17 ROI)
    double      neutronGunY_cm = -20.0;   // start position upstream of vessel tip

    // ── Gamma-source mode (event_type = 3, biased wall-background) ──
    // Re-emits capture-cascade gammas from a capture-vertex library CSV
    // (produced by scripts/make_capture_library.py from a neutron run).
    // CSV columns: volume,x_mm,y_mm,z_mm   (volume ∈ He3Cap_Al, He3Cap_CFRP)
    bool        gammaSourceMode = false;
    std::string captureLibFile;

    // ── Trajectory dump (event-display / cross-check) ───────
    // When set, writes per-step trajectories of neutrons and their charged
    // secondaries for the first trajDumpMaxEvents events to
    // <outFile>_traj[_t<tid>].csv.  Intended for neutron-path event displays
    // in the target; pair with a narrow --emin/--emax window for a
    // ~mono-energetic beam.  Cheap, but neutrons that fully thermalise take
    // thousands of elastic steps — keep the event count small.
    bool        trajDump          = false;
    int         trajDumpMaxEvents = 20;

    // ── Gas mixture ─────────────────────────────────────────
    std::string gas = "ArIso";

    // ── He-3 target dimensions (from STEP file: MASTINU X17 HPRV 00 01) ──
    // Gas bore: cylinder r=10 mm + r=10 mm hemispherical end caps (capsule shape)
    double he3_radius_cm      = 1.0;   // bore radius [cm]  (D=20 mm per STEP)
    double he3_half_length_cm = 2.0;   // half-length of cylinder section [cm]  (L=40 mm per STEP)

    // ── MM + trigger scint wall geometry ────────────────────
    // MM window front-face distance from origin.  Measured 2026-06-30: the
    // beam/target sit at the centre of the (roughly square) box formed by the
    // four inner mylar faces, so each face is half its opposing pair-span from
    // the origin — ±X arms (D,B): 40.8 cm span → 20.40 cm; ±Z arms (A,C):
    // 40.9 cm span → 20.45 cm.  (See GEOMETRY_COORDINATE_CONVENTION.md.)
    double mm_distance_x_cm = 20.40;  // ±X arm (D,B) window front face from origin [cm]
    double mm_distance_z_cm = 20.45;  // ±Z arm (A,C) window front face from origin [cm]

    // Per-MM tangential "pinwheel" shift [cm]: each MM is slid ⟂ its outward
    // normal in a clockwise (top-down) sense, along −uHat.  Indexed by sim arm:
    //   [0]=D(+X)  [1]=B(−X)  [2]=A(+Z)  [3]=C(−Z).
    // Measured 2026-06-30, halved 2026-07-14 (earlier values were 2× too large).
    double mm_pinwheel_shift_cm[4] = {1.55, 1.575, 1.635, 1.73};

    double mm_size_u_cm    = 38.0;   // MM active area: u [cm]
    double mm_size_v_cm    = 34.0;   // MM active area: v (along beam) [cm]

    // ══════════════════════════════════════════════════════════════════════
    //  Detector stack, inside → out (measured 2026-07-15; see
    //  GEOMETRY_CHANGE_CHECKLIST.md):  MM → SiPM wall → plastics → 1 LS layer.
    //  Depths below are quoted from the MM drift-mylar FRONT face (w = 0).
    //  The MM + PCB build is unchanged; the air gap after the PCB simply grows
    //  so the SiPM container front lands at the measured distance.
    // ══════════════════════════════════════════════════════════════════════

    // ── SiPM trigger-scintillator wall (50×50 cm, 20 bars of 2.5 cm) ────────
    // Centered on the mechanical STRUCTURE (u = 0), NOT on the pinwheel-shifted
    // MM.  Only 16 of the 20 bars are instrumented; the read-out window is
    // shifted 1 bar toward the MM (drop 3 bars on the far side, 1 on the MM
    // side).  Un-read bars are omitted from the sim (drawn transparent in the
    // diagrams).  Bars run along beam (v); array runs along the tangent (u).
    double sipm_front_from_mylar_cm = 11.0;  // mylar front → SiPM container front (measured)
    double sipm_container_depth_cm  = 3.5;   // SiPM container depth (measured 2026-07-17, all walls; was 3.3)
    double sipm_bar_width_cm        = 2.5;   // one bar width (u)
    int    sipm_n_bars              = 20;    // total bars across the 50 cm wall
    int    sipm_n_readout           = 16;    // instrumented bars (rest removed from sim)
    int    sipm_readout_shift_bars  = 1;     // read-out window shift toward the MM [bars]
    double sipm_scint_thick_cm      = 0.3;   // active scint depth (centered in container)
    double sipm_size_v_cm           = 50.0;  // bar length along beam (v)

    // ── Plastic scintillators (2 bars per arm, individually tape-wrapped) ───
    // Two 20×30 cm bars side-by-side in u, behind the SiPM wall.  Centered on
    // the MM (inherit the per-arm pinwheel shift) both horizontally (u) and
    // vertically (v = 0).
    double backscint_u_cm      = 20.0;  // each bar: u [cm]
    double backscint_v_cm      = 30.0;  // each bar: v (along beam) [cm]
    double backscint_thick_cm  = 2.5;   // each bar: depth [cm]  (measured ~2.5, nominal 2.0)
    double backscint_gap_cm    = 0.3;   // gap between wrapped bars [cm]
    double backscint_tape_um   = 200.0; // black mylar tape (outermost) [µm]
    double backscint_al_um     = 20.0;  // Al foil on scintillator surface [µm]
    // SiPM container back → plastics front, per arm [0]=D [1]=B [2]=A [3]=C
    // (measured 2026-07-17; was a uniform 7 cm placeholder)
    double gap_sipm_to_plastic_cm[4] = {6.5, 6.1, 6.3, 6.1};

    // ── Liquid scintillator vessel (1 per arm, STEP-derived) ───────────────
    // Source: STEP "LS X17.step" (Shapr3D export 2026-06-10, solid 'Corpo 03',
    // products SCINT1/SCINT2).  CFRP vessel: flat slab (45×45 cm face) +
    // 9 cm funnel lofting the slab cross-section down to a Ø50 mm neck where
    // the PMT is inserted halfway.  Liquid = slab + funnel + neck up to the
    // PMT face; filled with 6.5 L, so the slab faces bulge outward — modelled
    // as two ellipsoid domes whose height is solved from the fill volume
    // (≈13 mm/side).  STEP fine details dropped: r=3.5 mm fill tube along the
    // face, r=9 mm boss, r=4–10 mm edge fillets.  Wall = 2.6 mm CFRP
    // (2.0 structural + 0.6 liner; the 40 µm Al liner is neglected).
    // Placement surveyed 2026-07-17 (see the per-arm fields below):
    //   • Orientation: vessels VERTICAL on B and C (neck+PMT up, +v) and
    //     HORIZONTAL on A and D (neck+PMT along +u = "to the right when
    //     looking from behind the wall toward the target").
    //   • Depth: measured SiPM-container back → FLAT slab front face (taken
    //     at the vessel edge, off the bulge), per arm.
    //   • Height (v): from the bottom-bar survey — SiPM enclosures 62 cm
    //     tall with the SiPM bars assumed centred; a bar sits 6.8 cm above
    //     the enclosure bottom (DOUBLE-CHECK this 6.8); LS vessel bottoms
    //     measured above that bar: D 1.6, B 1.7, A 1.7, C 1.6 cm.  Chain:
    //     slab centre v = −31 + 6.8 + bottom + (slab half-height), where the
    //     half-height is 45.06/2 (vertical, B/C) or 45.12/2 (horizontal,
    //     A/D).  All four land within ±0.7 mm of v = 0.
    //   • Tangential (u), surveyed 2026-07-18: measured against a 60 cm
    //     horizontal reference centred on the SiPM wall (structure u = 0),
    //     left/right as seen from behind the wall (right = +u).  Left-edge
    //     distances (reference left edge → slab left edge) D/B/A/C =
    //     8.3/6.1/8.1/5.7 cm ⇒ slab centre u = d_L − 30 + half-width
    //     (22.53 horizontal A/D, 22.56 vertical B/C).  Right-edge
    //     cross-checks close the 60 cm to within 0.2–1.2 mm (B, C) and
    //     2.6–6.6 mm (D, A — likely the rounded slab corner); left trusted.
    //     B/C land on their pinwheel MM centres (≤2.4 mm); A/D sit ~2.3 cm
    //     right of theirs.
    double ls_slab_u_cm        = 45.12;  // slab outer width, u [cm]   (STEP 451.2 mm)
    double ls_slab_v_cm        = 45.06;  // slab outer length, v [cm]  (STEP 450.6 mm)
    double ls_slab_thick_cm    = 2.12;   // slab outer thickness, unbulged [cm] (STEP 21.2 mm)
    double ls_wall_mm          = 2.6;    // CFRP wall [mm] (2.0 structural + 0.6 liner)
    double ls_funnel_len_cm    = 9.0;    // funnel length: slab → neck [cm] (STEP 90 mm)
    double ls_neck_len_cm      = 12.97;  // neck cylinder length [cm]  (STEP 129.7 mm)
    double ls_neck_r_cm        = 2.5;    // neck outer radius [cm]     (STEP r=25 mm)
    double ls_fill_liters      = 6.5;    // liquid fill volume [L] → sets bulge dome height
    double ls_pmt_insert_frac  = 0.5;    // PMT window depth into the neck (fraction of neck)
    double ls_pmt_r_cm         = 2.2;    // PMT envelope radius [cm] (fits Ø44.8 neck bore)
    double ls_pmt_len_cm       = 11.5;   // PMT envelope length [cm]
    double ls_pmt_glass_mm     = 2.0;    // PMT borosilicate wall [mm]
    // Per-arm placement, arm order [0]=D [1]=B [2]=A [3]=C (measured 2026-07-17):
    // SiPM container back → FLAT slab front face (edge measurement, off the bulge):
    double ls_front_from_sipm_back_cm[4] = {12.3, 12.7, 12.3, 12.7};
    // Vessel rotation about the local w (depth) axis:
    //   0 = long axis vertical, neck/PMT up (+v);  −90 = horizontal, neck +u.
    double ls_rot_deg[4] = {-90.0, 0.0, -90.0, 0.0};
    // Slab-centre height relative to the beam (v = 0), from the bottom-bar
    // survey chain above:
    double ls_offset_v_cm[4] = {-0.04, +0.03, +0.06, -0.07};
    // Slab-centre u relative to the STRUCTURE (SiPM-wall centre, u = 0),
    // from the 60-cm-reference survey above (2026-07-18).  NB: the LS is
    // placed against the STRUCTURE, not the pinwheel-shifted MM.
    double ls_center_u_cm[4] = {+0.83, -1.34, +0.63, -1.74};
};
