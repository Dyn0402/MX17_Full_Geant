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
#include "G4DecayPhysics.hh"
#include "G4GenericBiasingPhysics.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList(double biasNCaptureFactor, double gammaCut_um)
    : G4VModularPhysicsList(), fGammaCut_um(gammaCut_um) {
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

    // Cross-section biasing for the rare ³He(n,γ) channel: wrap the neutron
    // nCapture process with a G4BiasingProcessInterface.  The actual biasing
    // (and its restriction to the He3Gas volume) is defined by the operator
    // attached in DetectorConstruction::ConstructSDandField().  Registered only
    // when biasing is requested so analog runs are byte-for-byte unchanged.
    if (biasNCaptureFactor > 1.0) {
        auto* biasing = new G4GenericBiasingPhysics();
        biasing->PhysicsBias("neutron", {"nCapture"});
        RegisterPhysics(biasing);
    }

    // NOTE: no G4NeutronTrackingCut. Its default 10 µs time limit kills slow
    // neutrons mid-flight (a 0.4 eV neutron needs ~20 µs to cross 20 cm),
    // silently removing nearly all sub-eV beam transport — fatal for the
    // thermal-capture physics this experiment depends on.  HP physics is
    // designed to track neutrons to thermalisation; let it.
}

void PhysicsList::SetCuts() {
    // Default range cuts -- 1 mm is G4 default
    // For gas detectors at low pressure we want to be more careful:
    // - Electrons: short range cut to capture delta rays in gas
    // - Photons  : default is fine
    SetCutValue(1.0 * mm,  "proton");
    SetCutValue(10.0 * um, "e-");       // short cut for electrons -- capture delta rays
    SetCutValue(10.0 * um, "e+");
    SetCutValue(fGammaCut_um * um, "gamma");  // default 100 um photon production threshold
}
