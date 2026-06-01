// ActionInitialization.cc

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "X17PrimaryGenerator.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(const SimConfig& cfg,
                                            const DetectorConstruction* detCon)
    : G4VUserActionInitialization(), fConfig(cfg), fDetCon(detCon) {}

void ActionInitialization::BuildForMaster() const {
    SetUserAction(new RunAction(fConfig, true));
}

void ActionInitialization::Build() const {
    auto* runAction   = new RunAction(fConfig, false);
    auto* eventAction = new EventAction(fConfig, runAction);

    SetUserAction(new X17PrimaryGenerator(fConfig, eventAction));
    SetUserAction(runAction);
    SetUserAction(eventAction);
    SetUserAction(new SteppingAction(fConfig, eventAction,
                                     fDetCon->GetArmAxes()));
}
