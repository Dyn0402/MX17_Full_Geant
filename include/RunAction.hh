#pragma once
// RunAction.hh
// Writes two ROOT TTrees per worker thread:
//   HitTree  — one entry per hit in a sensitive volume
//   EventTree — one entry per event (generated pair kinematics)
// Each thread writes its own file; merge with hadd.

#include "G4UserRunAction.hh"
#include "SimConfig.hh"
#include "EventData.hh"
#include <memory>
#include <string>

class G4Run;

class RunAction : public G4UserRunAction {
public:
    RunAction(const SimConfig& cfg, bool isMaster);
    ~RunAction() override;

    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction(const G4Run* run) override;

    void RecordEvent(const EventData& data);

private:
    const SimConfig& fConfig;
    bool             fIsMaster;

    struct Impl;
    std::unique_ptr<Impl> fImpl;

    long fTotalEvents = 0;
    long fTotalHits   = 0;
};
