#pragma once
// StackingAction.hh
// --roulette p1 p0 (NFS beam background): event-level Russian roulette.
// Stage 0 tracks the primary neutron alone; its secondaries wait.  When it is
// done, the arms the waiting secondaries can reach are predicted (straight
// line into a generous box per arm, or born inside one).  ≥ 2 arms: kept.
// 1 arm: kept with probability p1, 0 arms: p0, the survivors' event weight
// × 1/p.  Killed histories lose their stacks and are written with weight 0
// and no hits.  A whole-history roulette is unbiased for any observable,
// coincidences included, however crude the prediction is.

#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"
#include "SimConfig.hh"
#include <array>
#include <vector>

struct ArmAxes;

class StackingAction : public G4UserStackingAction {
public:
    StackingAction(const SimConfig& cfg, const std::array<ArmAxes, 4>& armAxes);
    ~StackingAction() override = default;

    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track) override;
    void NewStage() override;
    void PrepareNewEvent() override;

private:
    struct Sec { G4ThreeVector pos, dir; double ke; int pdg; };
    int ArmsReached(const Sec& s) const;   // bit mask of arms

    const SimConfig&              fConfig;
    const std::array<ArmAxes, 4>& fArmAxes;
    std::vector<Sec>              fSecs;
    int                           fStage = 0;
};
