// SteppingAction.cc
// Fires for every step. Records one HitData entry per step in a scored volume.
// The armID is the physical volume's copy number (set at placement time).
// Local u,v,w are computed from the arm's axis vectors stored in DetectorConstruction.

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "EventData.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4Threading.hh"

#include <cstring>
#include <sstream>

// Volumes scored — any step with edep > 0 in these gets a HitData entry
const std::map<std::string, bool> SteppingAction::kScoredVolumes = {
    {"DriftGas",    true},
    {"AmpGas",      true},
    {"PlasticScint",true},   // 3mm trigger scintillator
    {"LiqScint_1",  true},   // 2cm LAB layer 1
    {"LiqScint_2",  true},   // 2cm LAB layer 2
    {"BackScintL",  true},   // 2cm back plastic scint, left bar
    {"BackScintR",  true},   // 2cm back plastic scint, right bar
};

// Particles drawn by the trajectory dump: the neutron itself plus the charged
// products of (n,p)t and of any pair/cascade, so an event display shows the
// neutron path AND what it made at the end.  In neutron mode every track
// descends from the one primary neutron, so a particle-type filter suffices.
const std::set<std::string> SteppingAction::kTrajParticles = {
    "neutron", "proton", "triton", "alpha", "deuteron", "He3",
    "e-", "e+", "gamma",
};

SteppingAction::SteppingAction(const SimConfig& cfg, EventAction* eventAction,
                                const std::array<ArmAxes, 4>& armAxes)
    : G4UserSteppingAction(), fConfig(cfg),
      fEventAction(eventAction), fArmAxes(armAxes) {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    // ── Trajectory dump for event displays (--trajdump) ──────────────────
    if (fConfig.trajDump) DumpTrajectoryStep(step);

    // ── Neutron terminal-interaction recording (neutron mode) ────────────
    // Record where the PRIMARY neutron track ENDS via a hadronic interaction:
    // volume, position, and process name.  NOTE: ³He(n,p)t is filed under
    // "neutronInelastic" by the HP package — "nCapture" is radiative capture
    // only — so we must catch any process that kills the track.
    // Used for the capture budget and to build the gamma-source library
    // (library wants nCapture in Al/CFRP; budget wants everything).
    if (fConfig.neutronMode) {
        const G4Track* trk = step->GetTrack();
        if (trk->GetTrackID() == 1 &&
            trk->GetTrackStatus() == fStopAndKill &&
            trk->GetDefinition()->GetParticleName() == "neutron") {
            const G4StepPoint* post = step->GetPostStepPoint();
            const G4VProcess*  proc = post->GetProcessDefinedStep();
            const G4String procName = proc ? proc->GetProcessName() : "unknown";
            if (procName != "Transportation") {     // ignore world-boundary exit
                EventData& ed = fEventAction->GetEventData();
                const G4VPhysicalVolume* cpv = post->GetPhysicalVolume()
                    ? post->GetPhysicalVolume()
                    : step->GetPreStepPoint()->GetPhysicalVolume();
                if (cpv) ed.capture_vol = cpv->GetLogicalVolume()->GetName();
                ed.capture_proc = procName;
                G4ThreeVector cp = post->GetPosition();
                ed.cap_x = cp.x() / mm;
                ed.cap_y = cp.y() / mm;
                ed.cap_z = cp.z() / mm;
            }
        }
    }

    G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.0) return;

    const G4VPhysicalVolume* pv = step->GetPreStepPoint()->GetPhysicalVolume();
    if (!pv) return;

    // Volume name is the logical volume name (shared across arms)
    const std::string volName = pv->GetLogicalVolume()->GetName();
    if (kScoredVolumes.find(volName) == kScoredVolumes.end()) return;

    const G4Track* track = step->GetTrack();
    const G4StepPoint* pre = step->GetPreStepPoint();

    // For most volumes the arm copy number is on the physical volume itself.
    // BackScintL/R are nested: BackScint → BackScintAl → BackScintTape[arm copy]
    // Their own copy number is always 0 (wrong for arms 1-3), so we must walk
    // 2 levels up the touchable to reach the tape envelope that carries the
    // correct arm copy number.
    int armID;
    if (volName == "BackScintL" || volName == "BackScintR") {
        armID = pre->GetTouchable()->GetVolume(2)->GetCopyNo();
    } else {
        armID = pv->GetCopyNo();
    }
    if (armID < 0 || armID > 3) return;

    // Global hit position (midpoint of step)
    G4ThreeVector gpos = 0.5 * (pre->GetPosition() +
                                  step->GetPostStepPoint()->GetPosition());

    // Local u,v,w from arm axis vectors
    const ArmAxes& ax = fArmAxes[armID];
    G4ThreeVector delta = gpos - ax.frontFace;
    G4double u = delta.dot(ax.uHat) / mm;
    G4double v = delta.dot(ax.vHat) / mm;
    G4double w = delta.dot(ax.wHat) / mm;

    G4ThreeVector momDir = pre->GetMomentumDirection();

    HitData h;
    h.eventID  = fEventAction->GetEventData().eventID;
    h.trackID  = track->GetTrackID();
    h.parentID = track->GetParentID();
    h.armID    = armID;
    std::strncpy(h.detType,  volName.c_str(),                                      31); h.detType[31]  = '\0';
    std::strncpy(h.particle, track->GetDefinition()->GetParticleName().c_str(), 31); h.particle[31] = '\0';
    h.u    = u;
    h.v    = v;
    h.w    = w;
    h.edep = edep / eV;
    h.ke   = pre->GetKineticEnergy() / MeV;
    h.time = pre->GetGlobalTime()    / ns;
    h.gx   = gpos.x() / mm;
    h.gy   = gpos.y() / mm;
    h.gz   = gpos.z() / mm;
    h.px   = momDir.x();
    h.py   = momDir.y();
    h.pz   = momDir.z();

    fEventAction->GetEventData().hits.push_back(h);
}

// ─────────────────────────────────────────────────────────────────────────────
// One CSV row per step of a followed particle, for the first
// trajDumpMaxEvents events.  Segment endpoints (pre→post) are written
// explicitly so the consumer can draw each step as a line without stitching;
// KE at both ends shows the slowing-down along the path.
void SteppingAction::DumpTrajectoryStep(const G4Step* step) {
    const int eid = fEventAction->GetEventData().eventID;
    if (eid < 0 || eid >= fConfig.trajDumpMaxEvents) return;

    const G4Track* trk = step->GetTrack();
    const std::string pname = trk->GetDefinition()->GetParticleName();
    if (kTrajParticles.find(pname) == kTrajParticles.end()) return;

    if (!fTrajOpenTried) {
        fTrajOpenTried = true;
        std::ostringstream fn;
        fn << fConfig.outFile << "_traj";
        const int tid = G4Threading::G4GetThreadId();
        if (tid >= 0) fn << "_t" << tid;
        fn << ".csv";
        fTrajFile.open(fn.str());
        if (fTrajFile.is_open())
            fTrajFile << "eventID,trackID,parentID,particle,istep,"
                         "pre_x,pre_y,pre_z,post_x,post_y,post_z,"
                         "ke_pre_MeV,ke_post_MeV,volume,process,edep_MeV\n";
    }
    if (!fTrajFile.is_open()) return;

    const G4StepPoint* pre  = step->GetPreStepPoint();
    const G4StepPoint* post = step->GetPostStepPoint();
    const G4VProcess*  proc = post->GetProcessDefinedStep();
    const std::string  procName = proc ? proc->GetProcessName() : "none";
    const G4VPhysicalVolume* pv = pre->GetPhysicalVolume();
    const std::string  vol = pv ? pv->GetLogicalVolume()->GetName() : "OutOfWorld";
    const G4ThreeVector p0 = pre->GetPosition();
    const G4ThreeVector p1 = post->GetPosition();

    fTrajFile << eid << ',' << trk->GetTrackID() << ',' << trk->GetParentID() << ','
              << pname << ',' << trk->GetCurrentStepNumber() << ','
              << p0.x()/mm << ',' << p0.y()/mm << ',' << p0.z()/mm << ','
              << p1.x()/mm << ',' << p1.y()/mm << ',' << p1.z()/mm << ','
              << pre->GetKineticEnergy()/MeV << ',' << post->GetKineticEnergy()/MeV << ','
              << vol << ',' << procName << ','
              << step->GetTotalEnergyDeposit()/MeV << '\n';
}
