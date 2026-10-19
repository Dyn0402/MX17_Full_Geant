#pragma once
// BeamForceOperator.hh
// --force-beam (NFS beam background): forces the primary neutron's first
// collision inside the volume the generator picked for this event
// (EventTypeInfo::force_lv), with Geant4's truncated-exponential operation
// over all wrapped neutron processes: the collision point follows the analog
// e^(−Σx) law truncated at the volume exit, the process is chosen ∝ its cross
// section.  This is the "forced" branch of G4BOptrForceCollision (GB02)
// without the free-flight clone: the generator already sampled which volume
// is hit and carries the probability in the event weight, so one event stays
// one physical history.  After the collision (or once the primary leaves the
// volume) the operator steps aside and everything is analog.
//
// Attached to every logical volume except the World (DetectorConstruction).

#include "G4VBiasingOperator.hh"

class G4BOptnForceCommonTruncatedExp;
struct EventTypeInfo;

class BeamForceOperator : public G4VBiasingOperator {
public:
    BeamForceOperator();
    ~BeamForceOperator() override;

    void StartTracking(const G4Track* track) override;

private:
    G4VBiasingOperation* ProposeOccurenceBiasingOperation(
        const G4Track* track, const G4BiasingProcessInterface* callingProcess) override;
    G4VBiasingOperation* ProposeFinalStateBiasingOperation(
        const G4Track* track, const G4BiasingProcessInterface* callingProcess) override;
    G4VBiasingOperation* ProposeNonPhysicsBiasingOperation(
        const G4Track*, const G4BiasingProcessInterface*) override { return nullptr; }

    void OperationApplied(const G4BiasingProcessInterface* callingProcess,
                          G4BiasingAppliedCase biasingCase,
                          G4VBiasingOperation* operationApplied,
                          const G4VParticleChange* particleChangeProduced) override;
    void OperationApplied(const G4BiasingProcessInterface* callingProcess,
                          G4BiasingAppliedCase biasingCase,
                          G4VBiasingOperation* occurenceOperationApplied,
                          G4double weightForOccurenceInteraction,
                          G4VBiasingOperation* finalStateOperationApplied,
                          const G4VParticleChange* particleChangeProduced) override;
    void CheckOccured();

    G4BOptnForceCommonTruncatedExp* fForce;
    EventTypeInfo*                  fInfo = nullptr;   // current event, primary track only
    bool                            fInit = false;
};
