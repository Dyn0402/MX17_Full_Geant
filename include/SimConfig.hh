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

    // ── Gas mixture ─────────────────────────────────────────
    std::string gas = "ArIso";

    // ── He-3 target dimensions ──────────────────────────────
    double he3_radius_cm      = 1.5;   // gas cylinder radius [cm]  (diameter = 3 cm)
    double he3_half_length_cm = 4.0;   // half-length along Y [cm]  (total = 8 cm)

    // ── MM + trigger scint wall geometry ────────────────────
    double mm_distance_cm  = 22.0;   // MM window front face from origin [cm]
    double mm_size_u_cm    = 38.0;   // MM active area: u [cm]
    double mm_size_v_cm    = 34.0;   // MM active area: v (along beam) [cm]
    double scint_size_u_cm = 48.0;   // Trigger plastic scint: u [cm]
    double scint_size_v_cm = 48.0;   // Trigger plastic scint: v [cm]

    // ── Liquid scintillator stack (2 layers, 45×45 cm face) ─
    double ls_size_u_cm        = 45.0;  // LS active face: u [cm]
    double ls_size_v_cm        = 45.0;  // LS active face: v [cm]
    double ls_thick_cm         = 2.0;   // LAB layer thickness [cm]
    double ls_cfrp_mm          = 2.0;   // structural CFRP wall thickness [mm]
    // Inner liner before each LAB layer (inside the structural walls):
    double ls_inner_cfrp_um    = 600.0; // inner CFRP liner [µm]
    double ls_inner_al_um      = 40.0;  // Al liner [µm]

    // ── Back plastic scintillators (2 per arm, individually tape-wrapped) ──
    // Two 20×30 cm bars placed side-by-side in the u direction.
    double backscint_u_cm      = 20.0;  // each bar: u [cm]
    double backscint_v_cm      = 30.0;  // each bar: v (along beam) [cm]
    double backscint_thick_cm  = 2.0;   // each bar: depth [cm]
    double backscint_gap_cm    = 0.3;   // gap between wrapped bars [cm]
    double backscint_tape_um   = 200.0; // black mylar tape (outermost) [µm]
    double backscint_al_um     = 20.0;  // Al foil on scintillator surface [µm]

    // ── Clearances (air gaps) ──────────────────────────────
    double gap_pcb_to_scint_mm   = 20.0;  // PCB → trigger scint [mm]
    double gap_scint_to_ls_mm    = 20.0;  // trigger scint → LS stack [mm]
    double gap_ls_to_backscint_mm = 10.0; // LS stack → back plastic scints [mm]
};
