#pragma once
// EventData.hh — per-event data container shared between generator, stepping, and output

#include "HitData.hh"
#include "G4VUserEventInformation.hh"
#include <vector>

// Carried on the G4Event so generator → EventAction communication is MT-safe.
// (Raw pointers between user actions are fragile in Geant4 MT.)
struct EventTypeInfo : public G4VUserEventInformation {
    int    event_type    = 0;    // 0=X17, 1=IPC, -1=single
    double inv_mass_MeV  = 0.0;  // m_X17 or sampled Mee
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
    int              event_type = 0;  // 0 = X17, 1 = IPC, -1 = single particle
    PairKinematics   kin;
    std::vector<HitData> hits;

    void Reset() {
        eventID    = -1;
        event_type = 0;
        kin        = {};
        hits.clear();
    }
};
