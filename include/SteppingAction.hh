#pragma once
// SteppingAction.hh
// Records one HitData entry per step that deposits energy in a sensitive volume.

#include "G4UserSteppingAction.hh"
#include "SimConfig.hh"
#include <array>
#include <string>
#include <map>

class EventAction;
class G4Step;
struct ArmAxes;

class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction(const SimConfig& cfg, EventAction* eventAction,
                   const std::array<ArmAxes, 4>& armAxes);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step* step) override;

private:
    const SimConfig&                    fConfig;
    EventAction*                        fEventAction;
    const std::array<ArmAxes, 4>&       fArmAxes;

    // Volume names that are scored (detType string → true)
    static const std::map<std::string, bool> kScoredVolumes;
};
