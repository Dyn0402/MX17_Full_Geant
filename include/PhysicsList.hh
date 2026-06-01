#pragma once
// PhysicsList.hh
// Physics list for Micromegas simulation.
// Uses FTFP_BERT as hadronic base + G4EmStandardPhysics_option4 for EM
// (option4 = best EM, required for accurate low-energy electron/photon ionization)
// + G4RadioactiveDecay + G4StepLimiterPhysics

#include "G4VModularPhysicsList.hh"

class PhysicsList : public G4VModularPhysicsList {
public:
    PhysicsList();
    ~PhysicsList() override = default;
    void SetCuts() override;
};
