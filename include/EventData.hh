#pragma once
// EventData.hh — per-event data container shared between generator, stepping, and output

#include "HitData.hh"
#include "G4VUserEventInformation.hh"
#include <string>
#include <vector>

// Carried on the G4Event so generator → EventAction communication is MT-safe.
// (Raw pointers between user actions are fragile in Geant4 MT.)
struct EventTypeInfo : public G4VUserEventInformation {
    int    event_type    = 0;    // 0=X17, 1=IPC, -1=single, 2=neutron, 3=gamma-source
    double inv_mass_MeV  = 0.0;  // m_X17 or sampled Mee (modes 0/1); Eγ [MeV] (mode 3)
    double neutron_E_eV  = 0.0;  // primary neutron energy (mode 2)
    void Print() const override {}
};

struct PairKinematics {
    // Vertex (production point in He-3 gas) [mm]
    double vx = 0.0, vy = 0.0, vz = 0.0;
    // e- truth
    double em_ke  = 0.0;  // kinetic energy [MeV]
    double em_px  = 0.0, em_py = 0.0, em_pz = 0.0;  // momentum unit vector
    // e+ truth
    double ep_ke  = 0.0;
    double ep_px  = 0.0, ep_py = 0.0, ep_pz = 0.0;
    double openingAngle_deg = 0.0;
    double inv_mass_MeV = 0.0;  // pair invariant mass: m_X17 for signal, Mee for IPC [MeV]
};

struct EventData {
    int              eventID    = -1;
    int              event_type = 0;  // 0=X17, 1=IPC, -1=single, 2=neutron, 3=gamma-source
    PairKinematics   kin;
    std::vector<HitData> hits;

    // Neutron mode (event_type = 2): primary energy + terminal interaction.
    // capture_proc: "nCapture" = radiative; "neutronInelastic" = (n,p) etc.
    double      neutron_E_eV = 0.0;
    std::string capture_vol;                       // empty = escaped world
    std::string capture_proc;
    double      cap_x = 0.0, cap_y = 0.0, cap_z = 0.0;   // position [mm]

    void Reset() {
        eventID    = -1;
        event_type = 0;
        kin        = {};
        hits.clear();
        neutron_E_eV = 0.0;
        capture_vol.clear();
        capture_proc.clear();
        cap_x = cap_y = cap_z = 0.0;
    }
};
