#pragma once
// NCaptureBiasingOperator.hh
// Region-scoped cross-section biasing for the rare ³He(n,γ) radiative-capture
// channel (the X17/IPC source).  Attached to the He3Gas logical volume ONLY,
// it scales the nCapture interaction rate by a fixed factor so the channel is
// sampled ~factor× more often.  G4BOptnChangeCrossSection assigns each biased
// capture the weight 1/factor, so summing weights recovers the analog rate.
//
// Because the operator lives only on the gas volume, hydrogen/aluminium
// captures in the scintillators and capsule are never biased — neutron
// transport and the wall-background budget are untouched.  This is the
// standard GB06-style occurrence-biasing pattern.

#include "G4VBiasingOperator.hh"
#include "globals.hh"

class G4BOptnChangeCrossSection;
class G4BiasingProcessInterface;
class G4Track;

class NCaptureBiasingOperator : public G4VBiasingOperator {
public:
    explicit NCaptureBiasingOperator(G4double factor);
    ~NCaptureBiasingOperator() override;

private:
    // Occurrence biasing: scale the nCapture cross-section by fFactor.
    G4VBiasingOperation* ProposeOccurenceBiasingOperation(
        const G4Track* track,
        const G4BiasingProcessInterface* callingProcess) override;

    // We do not bias final states or non-physics processes.
    G4VBiasingOperation* ProposeFinalStateBiasingOperation(
        const G4Track*, const G4BiasingProcessInterface*) override { return nullptr; }
    G4VBiasingOperation* ProposeNonPhysicsBiasingOperation(
        const G4Track*, const G4BiasingProcessInterface*) override { return nullptr; }

    G4BOptnChangeCrossSection* fChangeCrossSection;
    G4double                   fFactor;
};
