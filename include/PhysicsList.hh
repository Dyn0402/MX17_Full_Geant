#pragma once
// PhysicsList.hh
// Physics list for Micromegas simulation.
// Uses FTFP_BERT as hadronic base + G4EmStandardPhysics_option4 for EM
// (option4 = best EM, required for accurate low-energy electron/photon ionization)
// + G4RadioactiveDecay + G4StepLimiterPhysics

#include "G4VModularPhysicsList.hh"

class PhysicsList : public G4VModularPhysicsList {
public:
    // biasNCaptureFactor > 1 wraps the neutron nCapture process for occurrence
    // biasing (the region-scoped operator is attached to the gas in
    // DetectorConstruction::ConstructSDandField).
    explicit PhysicsList(double biasNCaptureFactor = 1.0);
    ~PhysicsList() override = default;
    void SetCuts() override;
};
