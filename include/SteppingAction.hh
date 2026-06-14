#pragma once
// SteppingAction.hh
// Records one HitData entry per step that deposits energy in a sensitive volume.

#include "G4UserSteppingAction.hh"
#include "SimConfig.hh"
#include <array>
#include <fstream>
#include <set>
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
    // Trajectory dump (--trajdump): one CSV row per step of a neutron or one
    // of its charged secondaries, for the first trajDumpMaxEvents events.
    // File is per-thread (this action is constructed per worker), opened lazily.
    void DumpTrajectoryStep(const G4Step* step);

    const SimConfig&                    fConfig;
    EventAction*                        fEventAction;
    const std::array<ArmAxes, 4>&       fArmAxes;

    std::ofstream  fTrajFile;
    bool           fTrajOpenTried = false;

    // Volume names that are scored (detType string → true)
    static const std::map<std::string, bool> kScoredVolumes;
    // Particles followed by the trajectory dump
    static const std::set<std::string>       kTrajParticles;
};
