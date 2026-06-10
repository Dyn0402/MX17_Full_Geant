#pragma once
// X17PrimaryGenerator.hh
// Pair mode    : fires correlated e+/e- from X17→e+e- decay with correct Lorentz boost.
// Single mode  : fires one particle at a chosen energy/angle (efficiency cross-checks).
// Neutron mode : fires beam neutrons sampled from the EAR2 flux + radial profile.
// Gamma mode   : re-emits capture-cascade gammas from a capture-vertex library.

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
    G4double GenerateNeutron(G4Event* event);     // beam neutron, returns E_n [eV]
    G4double GenerateGammaSource(G4Event* event); // capture γ, returns Eγ [MeV]

    const SimConfig&          fConfig;
    std::unique_ptr<G4ParticleGun> fGun;

    G4ParticleDefinition* fElectron = nullptr;
    G4ParticleDefinition* fPositron = nullptr;
    G4ParticleDefinition* fNeutron  = nullptr;
    G4ParticleDefinition* fGamma    = nullptr;
};
