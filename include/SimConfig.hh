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
    double sipm_container_depth_cm  = 3.3;   // SiPM wall container depth; scint centered inside
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
    double gap_sipm_to_plastic_cm = 7.0; // SiPM container back → plastics front [cm] (MEASURE LATER)

    // ── Liquid scintillator (1 layer, 45×45 cm face) ───────────────────────
    // Single LAB layer inside a CFRP box (front wall + liner | LAB | rear wall).
    // Centered on the MM (assumed for now — MEASURE & UPDATE LATER).
    double ls_size_u_cm        = 45.0;  // LS active face: u [cm]
    double ls_size_v_cm        = 45.0;  // LS active face: v [cm]
    double ls_thick_cm         = 2.0;   // LAB layer thickness [cm]
    double ls_cfrp_mm          = 2.0;   // structural CFRP wall thickness [mm]
    double ls_inner_cfrp_um    = 600.0; // inner CFRP liner [µm]
    double ls_inner_al_um      = 40.0;  // Al liner [µm]
    double gap_plastic_to_ls_cm = 5.0;  // plastics back → liquid front [cm] (MEASURE LATER)
};
