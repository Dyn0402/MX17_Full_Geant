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
#include "EventData.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Exception.hh"
#include "Randomize.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <mutex>
#include <sstream>
#include <vector>

#ifdef USE_ROOT
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#endif

X17PrimaryGenerator::X17PrimaryGenerator(const SimConfig& cfg)
    : G4VUserPrimaryGeneratorAction(),
      fConfig(cfg), fGun(std::make_unique<G4ParticleGun>(1))
{
    G4ParticleTable* pt = G4ParticleTable::GetParticleTable();
    fElectron = pt->FindParticle("e-");
    fPositron = pt->FindParticle("e+");
    fNeutron  = pt->FindParticle("neutron");
    fGamma    = pt->FindParticle("gamma");
}

// ─────────────────────────────────────────────────────────────────────────────
void X17PrimaryGenerator::GeneratePrimaries(G4Event* event) {
    // EventTypeInfo is carried on the G4Event so EndOfEventAction can read it
    // without depending on a cross-action pointer (fragile in Geant4 MT).
    auto* info = new EventTypeInfo();

    if (fConfig.neutronMode) {
        info->event_type   = 2;
        info->neutron_E_eV = GenerateNeutron(event);
    } else if (fConfig.gammaSourceMode) {
        info->event_type   = 3;
        info->inv_mass_MeV = GenerateGammaSource(event);
    } else if (fConfig.singleParticle) {
        info->event_type = -1;
        GenerateSingle(event);
    } else if (G4UniformRand() < fConfig.ipc_fraction) {
        info->event_type   = 1;
        info->inv_mass_MeV = GenerateIPC(event);
    } else {
        info->event_type   = 0;
        info->inv_mass_MeV = fConfig.x17Mass_MeV;
        GeneratePair(event);
    }

    event->SetUserInformation(info);
}

// ─────────────────────────────────────────────────────────────────────────────
// Pair-vertex library (thermal self-shielding): He3Gas capture positions from
// a neutron run (make_capture_library.py --gas-lib).  When configured, pair
// vertices are drawn from these rows instead of uniformly in the gas.
namespace {
std::vector<std::array<double, 3>> gPairVtxLib;   // [mm]
std::once_flag gPairVtxOnce;

void LoadPairVertexLib(const SimConfig& cfg) {
    std::ifstream in(cfg.pairVertexLibFile);
    if (!in) {
        G4Exception("X17PrimaryGenerator", "PairVertexLib", FatalException,
                    ("Cannot open pair-vertex library: " + cfg.pairVertexLibFile).c_str());
        return;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.rfind("volume", 0) == 0)
            continue;
        std::istringstream ss(line);
        std::string vol, sx, sy, sz;
        if (!std::getline(ss, vol, ',') || !std::getline(ss, sx, ',') ||
            !std::getline(ss, sy, ',')  || !std::getline(ss, sz, ','))
            continue;
        if (vol.find("He3Gas") == std::string::npos) continue;
        gPairVtxLib.push_back({std::stod(sx), std::stod(sy), std::stod(sz)});
    }
    if (gPairVtxLib.empty())
        G4Exception("X17PrimaryGenerator", "PairVertexLib", FatalException,
                    "Pair-vertex library contains no He3Gas vertices");
    G4cout << "X17PrimaryGenerator: pair-vertex library — "
           << gPairVtxLib.size() << " He3Gas capture vertices" << G4endl;
}
}   // namespace

// Helper: sample a pair vertex.  Library mode (thermal absorption profile)
// when configured; otherwise uniform inside the He-3 gas capsule.
// Capsule = cylinder (r, ±Hl) + two hemispherical end caps of radius r.
// Uses rejection sampling from a bounding cylinder; acceptance ~89%.
static G4ThreeVector SampleHe3Vertex(const SimConfig& cfg) {
    if (!cfg.pairVertexLibFile.empty()) {
        std::call_once(gPairVtxOnce, LoadPairVertexLib, std::cref(cfg));
        const auto& p = gPairVtxLib[static_cast<size_t>(
            G4UniformRand() * gPairVtxLib.size()) % gPairVtxLib.size()];
        return G4ThreeVector(p[0] * mm, p[1] * mm, p[2] * mm);
    }
    const G4double R  = cfg.he3_radius_cm      * cm;
    const G4double Hl = cfg.he3_half_length_cm * cm;
    while (true) {
        G4double r    = R * std::sqrt(G4UniformRand());
        G4double phiV = CLHEP::twopi * G4UniformRand();
        G4double y    = (Hl + R) * (2.0 * G4UniformRand() - 1.0);
        G4double absY = std::abs(y);
        G4double dCap = absY - Hl;
        if (absY <= Hl || r*r + dCap*dCap <= R*R) {
            return G4ThreeVector(r * std::cos(phiV), y, r * std::sin(phiV));
        }
    }
}

// Helper: sample an isotropic unit direction.
static G4ThreeVector IsotropicDirection() {
    G4double cosT = 2.0 * G4UniformRand() - 1.0;
    G4double sinT = std::sqrt(1.0 - cosT * cosT);
    G4double phi  = CLHEP::twopi * G4UniformRand();
    return G4ThreeVector(sinT * std::cos(phi), cosT, sinT * std::sin(phi));
}

// ─────────────────────────────────────────────────────────────────────────────
// Core kinematics shared by both X17 and IPC.
// Fires e- (vertex 0) and e+ (vertex 1) from the given vertex.
// Truth is NOT stored in EventData here — EndOfEventAction reads it back from
// the G4Event's primary vertices, which is MT-safe.
static void GeneratePairFromParent(G4double m_parent, G4double E_parent,
                                   const G4ThreeVector& vertex,
                                   G4ParticleDefinition* electron,
                                   G4ParticleDefinition* positron,
                                   G4ParticleGun* gun,
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

    // ── Fire e- (vertex 0) then e+ (vertex 1) ────────────────────────────
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
    GeneratePairFromParent(m_x17, E_parent, vertex,
                           fElectron, fPositron, fGun.get(), event);
}

// ─────────────────────────────────────────────────────────────────────────────
// Returns the sampled invariant mass [MeV/c²] (stored in EventTypeInfo).
G4double X17PrimaryGenerator::GenerateIPC(G4Event* event) {
    const G4double me     = fElectron->GetPDGMass();
    G4double E_transition = fConfig.transition_energy_MeV * MeV;

    // Inverse CDF of dN/dMee ∝ 1/Mee on [2me, E_transition]
    G4double Mee = 2.0 * me * std::pow(E_transition / (2.0 * me), G4UniformRand());

    G4ThreeVector vertex = SampleHe3Vertex(fConfig);
    GeneratePairFromParent(Mee, E_transition, vertex,
                           fElectron, fPositron, fGun.get(), event);
    return Mee / MeV;
}

// ─────────────────────────────────────────────────────────────────────────────
// Neutron-beam mode (event_type = 2)
//
// Energy: sampled from `flux_n_pulse_NOisolet_100bpd` (neutrons/pulse per bin)
// of fluxEAR2-Ph3_in_different_units.root, restricted to
// [neutronEmin_eV, neutronEmax_eV]; log-uniform spread inside each bin.
//
// Transverse position: the radial profile row of Lambda2D
// (lamda2DvsEn_EAR2.root; 240 log E bins × 3000 r bins to 3 cm) nearest to the
// sampled energy.  Row weights are MC counts per radial bin — they already
// include the 2πr area factor, so the row CDF is sampled directly in r.
// Beam axis = +Y, so the transverse plane is world (x, z).
//
// The beam data are loaded once (master + worker threads share the static
// tables, guarded by std::once_flag) and used read-only afterwards.

#ifdef USE_ROOT
namespace {
struct BeamData {
    // energy CDF (per-bin lo/hi kept separately: window may contain
    // zero-flux bins, so accepted bins need not be contiguous)
    std::vector<double> eLo, eHi;  // accepted-bin edges [eV]
    std::vector<double> eCdf;      // cumulative flux, normalised to 1
    // radial profile per energy row
    std::vector<double> profLogE;  // log10(E/eV) row centres
    std::vector<std::vector<double>> rCdf;   // per-row CDF over radius
    std::vector<double> rEdges;    // radial bin edges [cm]
    bool ok = false;
};
BeamData      gBeam;
std::once_flag gBeamOnce;

void LoadBeamData(const SimConfig& cfg) {
    // ── flux ──────────────────────────────────────────────────────────────
    TFile ff(cfg.neutronFluxFile.c_str(), "READ");
    if (ff.IsZombie()) {
        G4Exception("X17PrimaryGenerator", "BeamData", FatalException,
                    ("Cannot open flux file: " + cfg.neutronFluxFile).c_str());
        return;
    }
    auto* hFlux = dynamic_cast<TH1*>(ff.Get("flux_n_pulse_NOisolet_100bpd"));
    if (!hFlux) {
        G4Exception("X17PrimaryGenerator", "BeamData", FatalException,
                    "Histogram flux_n_pulse_NOisolet_100bpd not found");
        return;
    }
    double cum = 0.0;
    for (int i = 1; i <= hFlux->GetNbinsX(); ++i) {
        double lo = hFlux->GetBinLowEdge(i);
        double hi = lo + hFlux->GetBinWidth(i);
        if (hi < cfg.neutronEmin_eV || lo > cfg.neutronEmax_eV) continue;
        double w = hFlux->GetBinContent(i);
        if (w <= 0) continue;
        gBeam.eLo.push_back(std::max(lo, cfg.neutronEmin_eV));
        gBeam.eHi.push_back(std::min(hi, cfg.neutronEmax_eV));
        cum += w;
        gBeam.eCdf.push_back(cum);
    }
    ff.Close();
    if (gBeam.eCdf.empty() || cum <= 0) {
        G4Exception("X17PrimaryGenerator", "BeamData", FatalException,
                    "No flux in the requested neutron energy window");
        return;
    }
    for (auto& c : gBeam.eCdf) c /= cum;

    // ── radial profile ────────────────────────────────────────────────────
    TFile fp(cfg.neutronProfileFile.c_str(), "READ");
    if (fp.IsZombie()) {
        G4Exception("X17PrimaryGenerator", "BeamData", FatalException,
                    ("Cannot open profile file: " + cfg.neutronProfileFile).c_str());
        return;
    }
    auto* h2 = dynamic_cast<TH2*>(fp.Get("Lambda2D"));
    if (!h2) {
        G4Exception("X17PrimaryGenerator", "BeamData", FatalException,
                    "Histogram Lambda2D not found");
        return;
    }
    const int nE = h2->GetNbinsX();
    const int nR = h2->GetNbinsY();
    const int rebin = 50;                       // 10 µm → 0.5 mm bins
    const int nRr = nR / rebin;
    gBeam.rEdges.resize(nRr + 1);
    for (int j = 0; j <= nRr; ++j)
        gBeam.rEdges[j] = h2->GetYaxis()->GetBinLowEdge(j * rebin + 1);

    for (int i = 1; i <= nE; ++i) {
        std::vector<double> cdf(nRr, 0.0);
        double rcum = 0.0;
        for (int j = 0; j < nRr; ++j) {
            double s = 0.0;
            for (int k = 1; k <= rebin; ++k)
                s += h2->GetBinContent(i, j * rebin + k);
            rcum += s;
            cdf[j] = rcum;
        }
        if (rcum <= 0) continue;                // empty energy row — skip
        for (auto& c : cdf) c /= rcum;
        gBeam.profLogE.push_back(std::log10(h2->GetXaxis()->GetBinCenter(i)));
        gBeam.rCdf.push_back(std::move(cdf));
    }
    fp.Close();
    if (gBeam.rCdf.empty()) {
        G4Exception("X17PrimaryGenerator", "BeamData", FatalException,
                    "Lambda2D has no populated energy rows");
        return;
    }
    gBeam.ok = true;
    G4cout << "X17PrimaryGenerator: beam data loaded — "
           << gBeam.eCdf.size() << " flux bins in ["
           << cfg.neutronEmin_eV << ", " << cfg.neutronEmax_eV << "] eV, "
           << gBeam.rCdf.size() << " profile rows" << G4endl;
}
}   // namespace
#endif  // USE_ROOT

G4double X17PrimaryGenerator::GenerateNeutron(G4Event* event) {
#ifndef USE_ROOT
    G4Exception("X17PrimaryGenerator", "NeutronMode", FatalException,
                "Neutron mode requires ROOT (USE_ROOT) to read the beam files");
    return 0.0;
#else
    std::call_once(gBeamOnce, LoadBeamData, std::cref(fConfig));

    // ── energy: inverse-CDF over flux bins, log-uniform inside the bin ───
    const double u = G4UniformRand();
    const auto it = std::lower_bound(gBeam.eCdf.begin(), gBeam.eCdf.end(), u);
    auto idx = static_cast<size_t>(it - gBeam.eCdf.begin());
    if (idx >= gBeam.eLo.size()) idx = gBeam.eLo.size() - 1;
    const double lo = gBeam.eLo[idx], hi = gBeam.eHi[idx];
    const double E_eV = lo * std::pow(hi / lo, G4UniformRand());

    // ── radius: nearest profile row in log10(E) ──────────────────────────
    const double le = std::log10(E_eV);
    const auto rit = std::lower_bound(gBeam.profLogE.begin(),
                                      gBeam.profLogE.end(), le);
    size_t row = static_cast<size_t>(rit - gBeam.profLogE.begin());
    if (row > 0 && (row == gBeam.profLogE.size() ||
                    le - gBeam.profLogE[row - 1] < gBeam.profLogE[row] - le))
        --row;
    const auto& cdf = gBeam.rCdf[row];
    const double ur = G4UniformRand();
    const auto jit = std::lower_bound(cdf.begin(), cdf.end(), ur);
    const auto j = static_cast<size_t>(jit - cdf.begin());
    const double r_cm = gBeam.rEdges[j] +
        G4UniformRand() * (gBeam.rEdges[j + 1] - gBeam.rEdges[j]);

    const double phi = CLHEP::twopi * G4UniformRand();
    const G4ThreeVector pos(r_cm * std::cos(phi) * cm,
                            fConfig.neutronGunY_cm * cm,
                            r_cm * std::sin(phi) * cm);

    fGun->SetParticleDefinition(fNeutron);
    fGun->SetParticlePosition(pos);
    fGun->SetParticleMomentumDirection(G4ThreeVector(0, 1, 0));
    fGun->SetParticleEnergy(E_eV * eV);
    fGun->GeneratePrimaryVertex(event);
    return E_eV;
#endif
}

// ─────────────────────────────────────────────────────────────────────────────
// Gamma-source mode (event_type = 3): biased wall-background generator.
// Re-emits one capture-cascade γ per event from a capture-vertex library
// (CSV: volume,x_mm,y_mm,z_mm — from scripts/make_capture_library.py).
// Line energies/intensities: IAEA-PGAA thermal-capture evaluations.
// Each event carries weight = (captures in volume per neutron) × (Σ I_γ)
//                              / N_generated — applied in analysis.

namespace {
struct CascadeLine { double e_MeV, intensity; };

// ²⁸Al prompt lines > 2 MeV (per capture); soft lines can't reach the trigger
const std::vector<CascadeLine> kAlLines = {
    {7.724, 0.275}, {7.693, 0.044}, {6.316, 0.019}, {6.102, 0.020},
    {5.134, 0.020}, {4.903, 0.019}, {4.734, 0.055}, {4.691, 0.026},
    {4.660, 0.029}, {4.259, 0.070}, {4.133, 0.068}, {3.876, 0.035},
    {3.591, 0.046}, {3.539, 0.072}, {3.465, 0.085}, {3.034, 0.082},
    {2.960, 0.036}, {2.821, 0.045}, {2.590, 0.035}, {2.283, 0.047},
};
// CFRP: ¹²C lines + ¹H 2.22 MeV weighted by capture share within CFRP
// (n_C·σ_C : n_H·σ_H ≈ 0.16 : 0.84 for 85/4 wt% C/H)
const std::vector<CascadeLine> kCfrpLines = {
    {4.945, 0.675 * 0.16}, {3.684, 0.323 * 0.16}, {2.223, 1.0 * 0.84},
};

struct CaptureLib {
    std::vector<std::array<double, 3>> posAl, posCfrp;   // [mm]
    bool ok = false;
};
CaptureLib    gCapLib;
std::once_flag gCapOnce;

void LoadCaptureLib(const SimConfig& cfg) {
    std::ifstream in(cfg.captureLibFile);
    if (!in) {
        G4Exception("X17PrimaryGenerator", "GammaSource", FatalException,
                    ("Cannot open capture library: " + cfg.captureLibFile).c_str());
        return;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.rfind("volume", 0) == 0)
            continue;
        std::istringstream ss(line);
        std::string vol, sx, sy, sz;
        if (!std::getline(ss, vol, ',') || !std::getline(ss, sx, ',') ||
            !std::getline(ss, sy, ',')  || !std::getline(ss, sz, ','))
            continue;
        std::array<double, 3> p{std::stod(sx), std::stod(sy), std::stod(sz)};
        if (vol.find("Al") != std::string::npos)        gCapLib.posAl.push_back(p);
        else if (vol.find("CFRP") != std::string::npos) gCapLib.posCfrp.push_back(p);
    }
    gCapLib.ok = !(gCapLib.posAl.empty() && gCapLib.posCfrp.empty());
    if (!gCapLib.ok)
        G4Exception("X17PrimaryGenerator", "GammaSource", FatalException,
                    "Capture library contains no Al/CFRP vertices");
    G4cout << "X17PrimaryGenerator: capture library — "
           << gCapLib.posAl.size() << " Al + "
           << gCapLib.posCfrp.size() << " CFRP vertices" << G4endl;
}

const CascadeLine& SampleLine(const std::vector<CascadeLine>& lines) {
    double tot = 0.0;
    for (const auto& l : lines) tot += l.intensity;
    double u = tot * G4UniformRand();
    for (const auto& l : lines) {
        u -= l.intensity;
        if (u <= 0) return l;
    }
    return lines.back();
}
}   // namespace

G4double X17PrimaryGenerator::GenerateGammaSource(G4Event* event) {
    std::call_once(gCapOnce, LoadCaptureLib, std::cref(fConfig));

    // pick a vertex proportional to library population (Al vs CFRP)
    const size_t nAl = gCapLib.posAl.size(), nCf = gCapLib.posCfrp.size();
    const bool useAl = G4UniformRand() * (nAl + nCf) < double(nAl);
    const auto& vec  = useAl ? gCapLib.posAl : gCapLib.posCfrp;
    const auto& p    = vec[static_cast<size_t>(G4UniformRand() * vec.size())
                           % vec.size()];
    const auto& line = SampleLine(useAl ? kAlLines : kCfrpLines);

    fGun->SetParticleDefinition(fGamma);
    fGun->SetParticlePosition(G4ThreeVector(p[0] * mm, p[1] * mm, p[2] * mm));
    fGun->SetParticleMomentumDirection(IsotropicDirection());
    fGun->SetParticleEnergy(line.e_MeV * MeV);
    fGun->GeneratePrimaryVertex(event);
    return line.e_MeV;
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
}
