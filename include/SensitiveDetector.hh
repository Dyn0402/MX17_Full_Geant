#pragma once
// SensitiveDetector.hh
// Minimal SD -- we do all scoring in SteppingAction for simplicity.
// This SD is required by G4 to register volumes as "active" but
// the actual hit collection is handled via SteppingAction + EventData.

#include "G4VSensitiveDetector.hh"
#include "SimConfig.hh"
#include <string>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class SensitiveDetector : public G4VSensitiveDetector {
public:
    SensitiveDetector(const std::string& name,
                      const std::string& hitsCollectionName,
                      const std::string& volumeLabel,
                      const SimConfig& cfg);
    ~SensitiveDetector() override = default;

    void   Initialize(G4HCofThisEvent* hce) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    void   EndOfEvent(G4HCofThisEvent* hce) override;

private:
    std::string      fVolumeLabel;
    const SimConfig& fConfig;
};
