#pragma once
// EventAction.hh

#include "G4UserEventAction.hh"
#include "SimConfig.hh"
#include "EventData.hh"

class RunAction;
class G4Event;

class EventAction : public G4UserEventAction {
public:
    EventAction(const SimConfig& cfg, RunAction* runAction);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    EventData& GetEventData() { return fData; }

private:
    const SimConfig& fConfig;
    RunAction*       fRunAction;
    EventData        fData;
};
