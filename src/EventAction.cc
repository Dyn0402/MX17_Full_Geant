// EventAction.cc

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(const SimConfig& cfg, RunAction* runAction)
    : G4UserEventAction(), fConfig(cfg), fRunAction(runAction) {}

void EventAction::BeginOfEventAction(const G4Event* event) {
    fData.Reset();
    fData.eventID = event->GetEventID();
}

void EventAction::EndOfEventAction(const G4Event*) {
    fRunAction->RecordEvent(fData);

    if (fConfig.verbose && (fData.eventID % 1000 == 0)) {
        G4cout << "[Event " << fData.eventID << "]"
               << "  nhits=" << fData.hits.size()
               << "  opening=" << fData.kin.openingAngle_deg << " deg"
               << G4endl;
    }
}
