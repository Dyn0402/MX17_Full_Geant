// BeamForceOperator.cc — see the header.
#include "BeamForceOperator.hh"
#include "EventData.hh"

#include "G4BOptnForceCommonTruncatedExp.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4BiasingProcessSharedData.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <cfloat>

BeamForceOperator::BeamForceOperator()
    : G4VBiasingOperator("BeamForce"),
      fForce(new G4BOptnForceCommonTruncatedExp("BeamForceTruncExp")) {}

BeamForceOperator::~BeamForceOperator() { delete fForce; }

void BeamForceOperator::StartTracking(const G4Track* track) {
    fInfo = nullptr;
    fInit = false;
    if (track->GetTrackID() != 1) return;
    const G4Event* ev = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    if (!ev) return;
    auto* info = static_cast<EventTypeInfo*>(ev->GetUserInformation());
    if (info && info->force_active) fInfo = info;
}

G4VBiasingOperation* BeamForceOperator::ProposeOccurenceBiasingOperation(
    const G4Track* track, const G4BiasingProcessInterface* callingProcess) {
    if (!fInfo || !fInfo->force_active) return nullptr;
    if (track->GetVolume()->GetLogicalVolume() != fInfo->force_lv) {
        // left the forced volume without colliding (should not happen): give up
        if (fInit) fInfo->force_active = false;
        return nullptr;
    }

    // As in G4BOptrForceCollision's "toBeForced" state: the first wrapped
    // process to be asked each step (re)initialises the shared operation.
    const G4bool first = callingProcess->GetIsFirstPostStepGPILInterface();
    if (first) {
        if (!fInit || fForce->GetInitialMomentum() != track->GetMomentum()) {
            fForce->Initialize(track);
            fInit = true;
        } else {
            fForce->UpdateForStep(track->GetStep());
        }
    }
    if (fForce->GetMaximumDistance() < DBL_MIN) {
        fInfo->force_active = false;
        return nullptr;
    }
    if (first) {
        const G4BiasingProcessSharedData* shared = callingProcess->GetSharedData();
        for (const auto* wrapper : shared->GetPhysicsBiasingProcessInterfaces()) {
            const G4double l = wrapper->GetWrappedProcess()->GetCurrentInteractionLength();
            if (l < DBL_MAX / 10.) fForce->AddCrossSection(wrapper->GetWrappedProcess(), 1.0 / l);
        }
        if (fForce->GetNumberOfSharing() > 0) fForce->Sample();
    }
    return callingProcess->GetWrappedProcess()->GetCurrentInteractionLength() < DBL_MAX / 10.
               ? fForce : nullptr;
}

G4VBiasingOperation* BeamForceOperator::ProposeFinalStateBiasingOperation(
    const G4Track*, const G4BiasingProcessInterface* callingProcess) {
    return callingProcess->GetCurrentOccurenceBiasingOperation();
}

void BeamForceOperator::CheckOccured() {
    if (fInfo && fInfo->force_active && fForce->GetInteractionOccured()) {
        fInfo->force_active = false;
        fInfo->force_done   = true;
    }
}

void BeamForceOperator::OperationApplied(const G4BiasingProcessInterface*, G4BiasingAppliedCase,
                                         G4VBiasingOperation*, const G4VParticleChange*) {
    CheckOccured();
}

void BeamForceOperator::OperationApplied(const G4BiasingProcessInterface*, G4BiasingAppliedCase,
                                         G4VBiasingOperation*, G4double,
                                         G4VBiasingOperation*, const G4VParticleChange*) {
    CheckOccured();
}
