#pragma once
// ActionInitialization.hh

#include "G4VUserActionInitialization.hh"
#include "SimConfig.hh"

class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization(const SimConfig& cfg,
                         const DetectorConstruction* detCon);
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    const SimConfig&            fConfig;
    const DetectorConstruction* fDetCon;
};
