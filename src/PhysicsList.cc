// PhysicsList.cc
// Builds the physics for the Micromegas simulation.
// Key considerations:
//   - Gammas  : photoelectric, Compton, pair production (all via EM option4)
//   - Neutrons: elastic + inelastic via FTFP_BERT (QGSP or FTFP + Bertini cascade)
//              + thermal neutrons via NeutronHP
//   - Electrons: accurate low-energy EM (Livermore model in option4 goes down to eV)
//   - Ionization: G4ionIonisation tracks energy loss in gas steps

#include "PhysicsList.hh"

// Modular physics components
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4DecayPhysics.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList() {
    SetVerboseLevel(0);

    // EM physics: option4 uses Livermore models below 100 keV for e-/gamma
    // This gives accurate photoelectric + Auger + ionization in low-Z gas
    RegisterPhysics(new G4EmStandardPhysics_option4(0));

    // EM extras: synchrotron, GDR etc (harmless to include)
    RegisterPhysics(new G4EmExtraPhysics(0));

    // Hadronic elastic with HP (high-precision neutron data, <20 MeV)
    RegisterPhysics(new G4HadronElasticPhysicsHP(0));

    // Hadronic inelastic with HP neutrons
    RegisterPhysics(new G4HadronPhysicsFTFP_BERT_HP(0));

    // Decay
    RegisterPhysics(new G4DecayPhysics(0));

    // Radioactive decay (useful for activation studies)
    RegisterPhysics(new G4RadioactiveDecayPhysics(0));

    // Step limiter (respects G4UserLimits set in detector volumes)
    RegisterPhysics(new G4StepLimiterPhysics());

    // Kill neutrons below 1 eV after tracking to avoid infinite loops
    RegisterPhysics(new G4NeutronTrackingCut(0));
}

void PhysicsList::SetCuts() {
    // Default range cuts -- 1 mm is G4 default
    // For gas detectors at low pressure we want to be more careful:
    // - Electrons: short range cut to capture delta rays in gas
    // - Photons  : default is fine
    SetCutValue(1.0 * mm,  "proton");
    SetCutValue(10.0 * um, "e-");       // short cut for electrons -- capture delta rays
    SetCutValue(10.0 * um, "e+");
    SetCutValue(0.1  * mm, "gamma");    // 0.1 mm photon production threshold
}
