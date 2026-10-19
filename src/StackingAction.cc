// StackingAction.cc — see the header.
#include "StackingAction.hh"
#include "DetectorConstruction.hh"
#include "EventData.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4StackManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "Randomize.hh"

#include <algorithm>
#include <bitset>
#include <cmath>

namespace {
// Arm box in the arm frame (front face centre, u across, v = beam, w outward):
// wider and deeper than any built stack (MM 400 × 360, big plastics 750 × 750,
// LS behind), so the prediction errs towards "reaches".
constexpr double kHalfU = 450.0, kHalfV = 450.0, kW0 = -10.0, kW1 = 700.0;   // mm

EventTypeInfo* CurrentInfo() {
    const G4Event* ev = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    return ev ? static_cast<EventTypeInfo*>(ev->GetUserInformation()) : nullptr;
}
}   // namespace

StackingAction::StackingAction(const SimConfig& cfg, const std::array<ArmAxes, 4>& armAxes)
    : G4UserStackingAction(), fConfig(cfg), fArmAxes(armAxes) {}

void StackingAction::PrepareNewEvent() {
    fSecs.clear();
    fStage = 0;
}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* track) {
    if (fStage > 0 || track->GetParentID() == 0) return fUrgent;
    fSecs.push_back({track->GetPosition() / mm, track->GetMomentumDirection(),
                     track->GetKineticEnergy() / MeV,
                     track->GetDefinition()->GetPDGEncoding()});
    return fWaiting;
}

int StackingAction::ArmsReached(const Sec& s) const {
    const int apdg = std::abs(s.pdg);
    // what can carry a signal out of where it was born
    bool travels = false;
    if (s.pdg == 22)                            travels = s.ke > 0.1;
    else if (apdg == 11)                        travels = s.ke > 0.5;
    else if (s.pdg == 2112)                     travels = true;
    else if (s.pdg == 2212 || s.pdg == 1000010020 || s.pdg == 1000010030)
                                                travels = s.ke > 2.0;
    int mask = 0;
    for (int a = 0; a < 4; ++a) {
        const ArmAxes& ax = fArmAxes[a];
        const G4ThreeVector d0 = s.pos - ax.frontFace / mm;
        const double p[3] = {d0.dot(ax.uHat), d0.dot(ax.vHat), d0.dot(ax.wHat)};
        const double lo[3] = {-kHalfU, -kHalfV, kW0}, hi[3] = {kHalfU, kHalfV, kW1};
        const bool inside = p[0] > lo[0] && p[0] < hi[0] && p[1] > lo[1] && p[1] < hi[1] &&
                            p[2] > lo[2] && p[2] < hi[2];
        if (inside) { mask |= 1 << a; continue; }
        if (!travels) continue;
        // slab test, ray p + t·d, t ≥ 0
        const double d[3] = {s.dir.dot(ax.uHat), s.dir.dot(ax.vHat), s.dir.dot(ax.wHat)};
        double t0 = 0.0, t1 = 1e30;
        bool hit = true;
        for (int k = 0; k < 3 && hit; ++k) {
            if (std::abs(d[k]) < 1e-12) { hit = p[k] > lo[k] && p[k] < hi[k]; continue; }
            double ta = (lo[k] - p[k]) / d[k], tb = (hi[k] - p[k]) / d[k];
            if (ta > tb) std::swap(ta, tb);
            t0 = std::max(t0, ta);
            t1 = std::min(t1, tb);
            hit = t0 < t1;
        }
        if (hit) mask |= 1 << a;
    }
    return mask;
}

void StackingAction::NewStage() {
    if (fStage > 0) return;
    fStage = 1;
    EventTypeInfo* info = CurrentInfo();
    int mask = 0;
    for (const auto& s : fSecs) mask |= ArmsReached(s);
    const int n = static_cast<int>(std::bitset<4>(mask).count());
    const double p = n >= 2 ? 1.0 : (n == 1 ? fConfig.rouletteP1 : fConfig.rouletteP0);
    if (info) info->n_arms_pred = n;
    if (p < 1.0 && G4UniformRand() >= p) {
        stackManager->clear();
        if (info) { info->killed = true; info->ev_weight = 0.0; }
        return;
    }
    if (info) info->ev_weight /= p;
    stackManager->ReClassify();
}
