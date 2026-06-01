// SensitiveDetector.cc
// The SD is intentionally lightweight -- scoring is in SteppingAction.

#include "SensitiveDetector.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"

SensitiveDetector::SensitiveDetector(const std::string& name,
                                     const std::string& hitsCollectionName,
                                     const std::string& volumeLabel,
                                     const SimConfig& cfg)
    : G4VSensitiveDetector(name),
      fVolumeLabel(volumeLabel),
      fConfig(cfg)
{
    collectionName.insert(hitsCollectionName);
}

void SensitiveDetector::Initialize(G4HCofThisEvent*) {}

G4bool SensitiveDetector::ProcessHits(G4Step*, G4TouchableHistory*) {
    // All hit processing is handled in SteppingAction for full cluster detail.
    return true;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*) {}
