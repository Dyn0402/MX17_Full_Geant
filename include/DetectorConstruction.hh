#pragma once
// DetectorConstruction.hh
// Builds: He-3 target at origin + 4 identical detector arms at ±X and ±Z.
// Beam axis = +Y (floor to ceiling).

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "SimConfig.hh"
#include <map>
#include <string>
#include <array>

class G4Material;

// Per-arm axis info for local-coordinate computation in SteppingAction
struct ArmAxes {
    G4ThreeVector frontFace; // world position of arm front face centre
    G4ThreeVector uHat;      // local u direction in world frame
    G4ThreeVector vHat;      // local v direction in world frame (always +Y)
    G4ThreeVector wHat;      // local w direction in world frame (depth, outward)
};

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    explicit DetectorConstruction(const SimConfig& cfg);
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Accessors for SteppingAction
    const std::array<ArmAxes, 4>& GetArmAxes() const { return fArmAxes; }

private:
    void DefineMaterials();
    G4Material* GetMat(const std::string& name) const;

    const SimConfig& fConfig;
    std::map<std::string, G4Material*> fMats;

    // Sensitive logical volumes (shared across all 4 arms via copy numbers)
    G4LogicalVolume* fDriftGasLV    = nullptr;
    G4LogicalVolume* fAmpGasLV      = nullptr;
    G4LogicalVolume* fPlScintLV     = nullptr;   // SiPM wall scint bar (placed 16×/arm)
    G4LogicalVolume* fLS1LV         = nullptr;   // single 2cm LAB layer
    // Back plastic scints: two bars per arm (L=left/−u, R=right/+u in arm frame)
    // Each placed 4× (copyNo = armID); separate LVs so volume name encodes L/R.
    G4LogicalVolume* fBackScintLLV  = nullptr;
    G4LogicalVolume* fBackScintRLV  = nullptr;

    G4LogicalVolume* fHe3GasLV      = nullptr;   // target gas: nCapture-bias host

    std::array<ArmAxes, 4> fArmAxes;
};
