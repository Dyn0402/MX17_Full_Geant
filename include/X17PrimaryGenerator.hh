#pragma once
// X17PrimaryGenerator.hh
// Pair mode   : fires correlated e+/e- from X17→e+e- decay with correct Lorentz boost.
// Single mode : fires one particle at a chosen energy/angle (efficiency cross-checks).

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "SimConfig.hh"
#include "EventData.hh"
#include <memory>

class G4Event;
class G4ParticleDefinition;

class X17PrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
    explicit X17PrimaryGenerator(const SimConfig& cfg);
    ~X17PrimaryGenerator() override = default;

    void GeneratePrimaries(G4Event* event) override;

private:
    void     GeneratePair(G4Event* event);   // X17 → e+e-
    G4double GenerateIPC(G4Event* event);    // IPC γ* → e+e-, returns Mee [MeV]
    void     GenerateSingle(G4Event* event); // single particle (cross-check mode)

    const SimConfig&          fConfig;
    std::unique_ptr<G4ParticleGun> fGun;

    G4ParticleDefinition* fElectron = nullptr;
    G4ParticleDefinition* fPositron = nullptr;
};
