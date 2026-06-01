// X17PrimaryGenerator.cc
//
// Reaction: 3He + n  →  4He*  (at rest in the lab, slow neutron capture)
// De-excitation of 4He* at transition_energy_MeV (~20.58 MeV), two modes:
//
//   X17 mode  (event_type = 0):
//     4He* → 4He_gs + X17  (m_X17 = 16.8 MeV, emitted isotropically)
//     X17 → e+ e-  (isotropic decay in X17 rest frame, boosted to lab)
//
//   IPC mode  (event_type = 1):
//     4He* → 4He_gs + γ*  (virtual photon, emitted isotropically)
//     Invariant mass Mee sampled from dN/dMee ∝ 1/Mee  (log-uniform, 2me → E_transition)
//     γ* → e+ e-  (isotropic in γ* rest frame, boosted to lab)
//     Same code structure as X17 — only the "parent mass" differs.
//
//   Single mode  (event_type = -1):
//     One particle at fixed energy/direction (efficiency cross-checks).
//
// The boost direction is always randomly isotropic.  The boost magnitude is:
//   β = sqrt(E_parent² - m_parent²) / E_parent
// where E_parent ≈ transition_energy (nuclear recoil is negligible for 4He).

#include "X17PrimaryGenerator.hh"
#include "EventAction.hh"
#include "EventData.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include <cmath>

X17PrimaryGenerator::X17PrimaryGenerator(const SimConfig& cfg,
                                          EventAction* eventAction)
    : G4VUserPrimaryGeneratorAction(),
      fConfig(cfg), fEventAction(eventAction),
      fGun(std::make_unique<G4ParticleGun>(1))
{
    G4ParticleTable* pt = G4ParticleTable::GetParticleTable();
    fElectron = pt->FindParticle("e-");
    fPositron = pt->FindParticle("e+");
}

// ─────────────────────────────────────────────────────────────────────────────
void X17PrimaryGenerator::GeneratePrimaries(G4Event* event) {
    EventData& ed = fEventAction->GetEventData();
    if (fConfig.singleParticle) {
        ed.event_type = -1;
        GenerateSingle(event);
    } else if (G4UniformRand() < fConfig.ipc_fraction) {
        ed.event_type = 1;
        GenerateIPC(event);
    } else {
        ed.event_type = 0;
        GeneratePair(event);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: sample a uniform random vertex inside the He-3 cylinder.
static G4ThreeVector SampleHe3Vertex(const SimConfig& cfg) {
    G4double R  = cfg.he3_radius_cm      * cm;
    G4double Hl = cfg.he3_half_length_cm * cm;
    G4double r    = R * std::sqrt(G4UniformRand());
    G4double phiV = CLHEP::twopi * G4UniformRand();
    return G4ThreeVector(r * std::cos(phiV),
                         Hl * (2.0 * G4UniformRand() - 1.0),
                         r * std::sin(phiV));
}

// Helper: sample an isotropic unit direction.
static G4ThreeVector IsotropicDirection() {
    G4double cosT = 2.0 * G4UniformRand() - 1.0;
    G4double sinT = std::sqrt(1.0 - cosT * cosT);
    G4double phi  = CLHEP::twopi * G4UniformRand();
    return G4ThreeVector(sinT * std::cos(phi), cosT, sinT * std::sin(phi));
}

// ─────────────────────────────────────────────────────────────────────────────
// Core kinematics shared by both X17 and IPC:
//   parent mass m_parent, parent energy E_parent (both in Geant4 internal units).
// Fills EventData, fires both particles from the given vertex.
static void GeneratePairFromParent(G4double m_parent, G4double E_parent,
                                   const G4ThreeVector& vertex,
                                   G4ParticleDefinition* electron,
                                   G4ParticleDefinition* positron,
                                   G4ParticleGun* gun,
                                   EventData& ed,
                                   G4Event* event)
{
    const G4double me = electron->GetPDGMass();

    // ── Boost vector: parent emitted isotropically from 4He* at rest ─────
    G4double p_parent = std::sqrt(E_parent * E_parent - m_parent * m_parent);
    G4ThreeVector boostVec = (p_parent / E_parent) * IsotropicDirection();

    // ── Decay in parent rest frame: e+e- back-to-back isotropically ──────
    G4double E_e = m_parent / 2.0;
    G4double p_e = (E_e > me) ? std::sqrt(E_e * E_e - me * me) : 0.0;
    G4ThreeVector decDir = IsotropicDirection();

    G4LorentzVector p4em(-p_e * decDir, E_e);
    G4LorentzVector p4ep( p_e * decDir, E_e);

    // ── Boost to lab ──────────────────────────────────────────────────────
    p4em.boost(boostVec);
    p4ep.boost(boostVec);

    G4ThreeVector momDir_em = p4em.vect().unit();
    G4ThreeVector momDir_ep = p4ep.vect().unit();
    G4double ke_em = std::max(0.0, p4em.e() - me);
    G4double ke_ep = std::max(0.0, p4ep.e() - me);

    G4double cosOpen = momDir_em.dot(momDir_ep);
    G4double openDeg = std::acos(std::max(-1.0, std::min(1.0, cosOpen))) / deg;

    // ── Store truth ───────────────────────────────────────────────────────
    ed.kin.vx  = vertex.x() / mm;
    ed.kin.vy  = vertex.y() / mm;
    ed.kin.vz  = vertex.z() / mm;
    ed.kin.em_ke = ke_em / MeV;
    ed.kin.ep_ke = ke_ep / MeV;
    ed.kin.em_px = momDir_em.x(); ed.kin.em_py = momDir_em.y(); ed.kin.em_pz = momDir_em.z();
    ed.kin.ep_px = momDir_ep.x(); ed.kin.ep_py = momDir_ep.y(); ed.kin.ep_pz = momDir_ep.z();
    ed.kin.openingAngle_deg = openDeg;
    ed.kin.inv_mass_MeV     = m_parent / MeV;

    // ── Fire particles ────────────────────────────────────────────────────
    gun->SetParticlePosition(vertex);

    gun->SetParticleDefinition(electron);
    gun->SetParticleEnergy(ke_em);
    gun->SetParticleMomentumDirection(momDir_em);
    gun->GeneratePrimaryVertex(event);

    gun->SetParticleDefinition(positron);
    gun->SetParticleEnergy(ke_ep);
    gun->SetParticleMomentumDirection(momDir_ep);
    gun->GeneratePrimaryVertex(event);
}

// ─────────────────────────────────────────────────────────────────────────────
void X17PrimaryGenerator::GeneratePair(G4Event* event) {
    G4double m_x17    = fConfig.x17Mass_MeV          * MeV;
    G4double E_parent = fConfig.transition_energy_MeV * MeV;
    G4ThreeVector vertex = SampleHe3Vertex(fConfig);
    EventData& ed = fEventAction->GetEventData();
    GeneratePairFromParent(m_x17, E_parent, vertex,
                           fElectron, fPositron, fGun.get(), ed, event);
}

// ─────────────────────────────────────────────────────────────────────────────
void X17PrimaryGenerator::GenerateIPC(G4Event* event) {
    const G4double me     = fElectron->GetPDGMass();
    G4double E_transition = fConfig.transition_energy_MeV * MeV;

    // ── Sample virtual photon invariant mass from dN/dMee ∝ 1/Mee ────────
    // Inverse CDF of 1/Mee on [2me, E_transition]: Mee = 2me * (E_t/2me)^U
    G4double Mee = 2.0 * me * std::pow(E_transition / (2.0 * me), G4UniformRand());

    G4ThreeVector vertex = SampleHe3Vertex(fConfig);
    EventData& ed = fEventAction->GetEventData();
    GeneratePairFromParent(Mee, E_transition, vertex,
                           fElectron, fPositron, fGun.get(), ed, event);
}

// ─────────────────────────────────────────────────────────────────────────────
void X17PrimaryGenerator::GenerateSingle(G4Event* event) {
    G4ParticleTable* pt = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pd = pt->FindParticle(fConfig.singleParticleName);
    if (!pd) pd = fElectron;

    G4double theta = fConfig.singleParticleTheta_deg * deg;
    G4double phi   = fConfig.singleParticlePhi_deg   * deg;
    G4ThreeVector momDir(std::sin(theta) * std::cos(phi),
                          std::cos(theta),
                          std::sin(theta) * std::sin(phi));

    fGun->SetParticlePosition(G4ThreeVector(0, 0, 0));
    fGun->SetParticleDefinition(pd);
    fGun->SetParticleEnergy(fConfig.singleParticleEnergy_MeV * MeV);
    fGun->SetParticleMomentumDirection(momDir);
    fGun->GeneratePrimaryVertex(event);

    EventData& ed = fEventAction->GetEventData();
    ed.kin.vx = ed.kin.vy = ed.kin.vz = 0.0;
    ed.kin.em_ke = fConfig.singleParticleEnergy_MeV;
    ed.kin.em_px = momDir.x(); ed.kin.em_py = momDir.y(); ed.kin.em_pz = momDir.z();
    ed.kin.inv_mass_MeV = 0.0;
}
