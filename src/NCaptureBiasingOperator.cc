// NCaptureBiasingOperator.cc
#include "NCaptureBiasingOperator.hh"

#include "G4BOptnChangeCrossSection.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4VProcess.hh"
#include "G4Track.hh"

#include <cfloat>

NCaptureBiasingOperator::NCaptureBiasingOperator(G4double factor)
    : G4VBiasingOperator("NCaptureXSBias"),
      fChangeCrossSection(nullptr),
      fFactor(factor) {
    fChangeCrossSection = new G4BOptnChangeCrossSection("nCaptureXSChange");
}

NCaptureBiasingOperator::~NCaptureBiasingOperator() {
    delete fChangeCrossSection;
}

G4VBiasingOperation* NCaptureBiasingOperator::ProposeOccurenceBiasingOperation(
    const G4Track* /*track*/,
    const G4BiasingProcessInterface* callingProcess) {

    // Only the radiative-capture process is boosted (belt-and-suspenders: only
    // nCapture is wrapped for biasing in the physics list anyway).
    if (callingProcess->GetWrappedProcess()->GetProcessName() != "nCapture")
        return nullptr;

    // The wrapped process has already computed its analog interaction length at
    // the current point.  A near-infinite length means "no interaction here" —
    // leave it analog.
    G4double analogInteractionLength =
        callingProcess->GetWrappedProcess()->GetCurrentInteractionLength();
    if (analogInteractionLength > DBL_MAX / 10.0) return nullptr;

    G4double analogXS = 1.0 / analogInteractionLength;
    G4double biasedXS = fFactor * analogXS;

    fChangeCrossSection->SetBiasedCrossSection(biasedXS);
    fChangeCrossSection->Sample();
    return fChangeCrossSection;
}
