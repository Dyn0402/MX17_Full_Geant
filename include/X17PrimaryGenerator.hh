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
class EventAction;

class X17PrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
    X17PrimaryGenerator(const SimConfig& cfg, EventAction* eventAction);
    ~X17PrimaryGenerator() override = default;

    void GeneratePrimaries(G4Event* event) override;

private:
    void GeneratePair(G4Event* event);    // X17 → e+e-
    void GenerateIPC(G4Event* event);     // IPC γ* → e+e- (virtual photon, sampled Mee)
    void GenerateSingle(G4Event* event);  // single particle (cross-check mode)

    const SimConfig&          fConfig;
    EventAction*              fEventAction;
    std::unique_ptr<G4ParticleGun> fGun;

    G4ParticleDefinition* fElectron = nullptr;
    G4ParticleDefinition* fPositron = nullptr;
};
