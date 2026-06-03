// EventAction.cc

#include "EventAction.hh"
#include "RunAction.hh"
#include "EventData.hh"

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LorentzVector.hh"

#include <cmath>
#include <algorithm>

EventAction::EventAction(const SimConfig& cfg, RunAction* runAction)
    : G4UserEventAction(), fConfig(cfg), fRunAction(runAction) {}

void EventAction::BeginOfEventAction(const G4Event* event) {
    fData.Reset();
    fData.eventID = event->GetEventID();
}

void EventAction::EndOfEventAction(const G4Event* event) {
    // ── Read event_type from EventTypeInfo carried on the G4Event ─────────
    // Using the G4Event as the communication channel is MT-safe: both
    // GeneratePrimaries and EndOfEventAction receive the same G4Event*.
    auto* info = static_cast<EventTypeInfo*>(event->GetUserInformation());
    if (info) {
        fData.event_type = info->event_type;
    }

    // ── Read truth kinematics from primary vertices ───────────────────────
    // G4ParticleGun::GeneratePrimaryVertex adds one vertex per call.
    // Pair events: vertex 0 = e-, vertex 1 = e+.
    G4int nVtx = event->GetNumberOfPrimaryVertex();
    if (nVtx >= 2) {
        G4PrimaryVertex*   v0 = event->GetPrimaryVertex(0);
        G4PrimaryVertex*   v1 = event->GetPrimaryVertex(1);
        G4PrimaryParticle* p0 = v0->GetPrimary();
        G4PrimaryParticle* p1 = v1->GetPrimary();

        const G4double me = electron_mass_c2;  // 0.511 MeV (CLHEP)

        fData.kin.vx = v0->GetX0() / mm;
        fData.kin.vy = v0->GetY0() / mm;
        fData.kin.vz = v0->GetZ0() / mm;

        G4ThreeVector mom0 = p0->GetMomentum();
        G4ThreeVector mom1 = p1->GetMomentum();
        G4double E0 = std::sqrt(mom0.mag2() + me * me);
        G4double E1 = std::sqrt(mom1.mag2() + me * me);

        fData.kin.em_ke = (E0 - me) / MeV;
        fData.kin.ep_ke = (E1 - me) / MeV;

        G4ThreeVector d0 = mom0.unit();
        G4ThreeVector d1 = mom1.unit();
        fData.kin.em_px = d0.x(); fData.kin.em_py = d0.y(); fData.kin.em_pz = d0.z();
        fData.kin.ep_px = d1.x(); fData.kin.ep_py = d1.y(); fData.kin.ep_pz = d1.z();

        G4double cosOpen = d0.dot(d1);
        fData.kin.openingAngle_deg =
            std::acos(std::max(-1.0, std::min(1.0, cosOpen))) / deg;

        // Invariant mass from 4-momenta
        G4LorentzVector lv0(mom0, E0), lv1(mom1, E1);
        fData.kin.inv_mass_MeV = (lv0 + lv1).mag() / MeV;

    } else if (nVtx == 1) {
        // Single-particle mode
        G4PrimaryVertex*   v0 = event->GetPrimaryVertex(0);
        G4PrimaryParticle* p0 = v0->GetPrimary();
        const G4double me = electron_mass_c2;
        G4double p0mag = p0->GetMomentum().mag();

        fData.kin.vx    = v0->GetX0() / mm;
        fData.kin.vy    = v0->GetY0() / mm;
        fData.kin.vz    = v0->GetZ0() / mm;
        fData.kin.em_ke = (std::sqrt(p0mag * p0mag + me * me) - me) / MeV;
        G4ThreeVector d0 = p0->GetMomentum().unit();
        fData.kin.em_px = d0.x(); fData.kin.em_py = d0.y(); fData.kin.em_pz = d0.z();
    }

    fRunAction->RecordEvent(fData);

    if (fConfig.verbose && (fData.eventID % 1000 == 0)) {
        G4cout << "[Event " << fData.eventID << "]"
               << "  type=" << fData.event_type
               << "  nhits=" << fData.hits.size()
               << "  opening=" << fData.kin.openingAngle_deg << " deg"
               << "  em_ke=" << fData.kin.em_ke << " MeV"
               << G4endl;
    }
}
