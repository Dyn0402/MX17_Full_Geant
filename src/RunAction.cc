// RunAction.cc
// Per-worker-thread ROOT output:
//   HitTree  — one entry per hit in a sensitive volume
//   EventTree — one entry per event (generated pair kinematics)
// Merge files with: hadd merged.root x17_output_t*.root

#include "RunAction.hh"
#include "G4Run.hh"
#include "G4Threading.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdio>

#ifdef USE_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

struct RunAction::Impl {
#ifdef USE_ROOT
    TFile* rootFile  = nullptr;
    TTree* hitTree   = nullptr;
    TTree* evtTree   = nullptr;

    // HitTree branches
    Int_t    h_eventID, h_trackID, h_parentID, h_armID;
    Char_t   h_detType[32], h_particle[32];
    Double_t h_u, h_v, h_w;
    Double_t h_edep, h_ke, h_time;
    Double_t h_gx, h_gy, h_gz;
    Double_t h_px, h_py, h_pz;

    // EventTree branches
    Int_t    e_eventID;
    Int_t    e_event_type;                  // 0=X17, 1=IPC, -1=single, 2=neutron, 3=gamma
    Double_t e_vtx_x, e_vtx_y, e_vtx_z;   // production vertex [mm]
    Double_t e_inv_mass;                    // pair invariant mass [MeV]: m_X17 or Mee; Eγ (mode 3)
    Double_t e_em_ke, e_ep_ke;
    Double_t e_em_px, e_em_py, e_em_pz;
    Double_t e_ep_px, e_ep_py, e_ep_pz;
    Double_t e_openingAngle;
    // Neutron mode (event_type = 2)
    Double_t e_neutron_E_eV;                // primary neutron energy [eV]
    Char_t   e_capture_vol[32];             // logical volume of first capture ("" = none)
    Double_t e_cap_x, e_cap_y, e_cap_z;     // capture position [mm]
#else
    std::ofstream hitFile;
    std::ofstream evtFile;
#endif
};

RunAction::RunAction(const SimConfig& cfg, bool isMaster)
    : G4UserRunAction(), fConfig(cfg), fIsMaster(isMaster),
      fImpl(std::make_unique<Impl>()) {}

RunAction::~RunAction() = default;

// ─────────────────────────────────────────────────────────────────────────────
void RunAction::BeginOfRunAction(const G4Run*) {
    fTotalEvents = fTotalHits = 0;
    if (fIsMaster) return;

    std::ostringstream ss;
    ss << fConfig.outFile;
    int tid = G4Threading::G4GetThreadId();
    if (tid >= 0) ss << "_t" << tid;

#ifdef USE_ROOT
    std::string fname = ss.str() + ".root";
    // Delete any existing file first — ROOT's RECREATE may not truncate properly
    // over EOS FUSE or other networked filesystems, causing stale cycles to persist.
    std::remove(fname.c_str());
    fImpl->rootFile = TFile::Open(fname.c_str(), "RECREATE");
    if (!fImpl->rootFile || fImpl->rootFile->IsZombie()) {
        G4cerr << "ERROR: Cannot open ROOT file " << fname << G4endl;
        return;
    }

    // ── HitTree ───────────────────────────────────────────
    fImpl->hitTree = new TTree("HitTree", "Per-hit data");
    fImpl->hitTree->Branch("eventID",  &fImpl->h_eventID);
    fImpl->hitTree->Branch("trackID",  &fImpl->h_trackID);
    fImpl->hitTree->Branch("parentID", &fImpl->h_parentID);
    fImpl->hitTree->Branch("armID",    &fImpl->h_armID);
    fImpl->hitTree->Branch("detType",  fImpl->h_detType,  "detType[32]/C");
    fImpl->hitTree->Branch("particle", fImpl->h_particle, "particle[32]/C");
    fImpl->hitTree->Branch("u",        &fImpl->h_u);
    fImpl->hitTree->Branch("v",        &fImpl->h_v);
    fImpl->hitTree->Branch("w",        &fImpl->h_w);
    fImpl->hitTree->Branch("edep",     &fImpl->h_edep);
    fImpl->hitTree->Branch("ke",       &fImpl->h_ke);
    fImpl->hitTree->Branch("time",     &fImpl->h_time);
    fImpl->hitTree->Branch("gx",       &fImpl->h_gx);
    fImpl->hitTree->Branch("gy",       &fImpl->h_gy);
    fImpl->hitTree->Branch("gz",       &fImpl->h_gz);
    fImpl->hitTree->Branch("px",       &fImpl->h_px);
    fImpl->hitTree->Branch("py",       &fImpl->h_py);
    fImpl->hitTree->Branch("pz",       &fImpl->h_pz);

    // ── EventTree ─────────────────────────────────────────
    fImpl->evtTree = new TTree("EventTree", "Per-event generated kinematics");
    fImpl->evtTree->Branch("eventID",      &fImpl->e_eventID);
    fImpl->evtTree->Branch("event_type",   &fImpl->e_event_type); // 0=X17, 1=IPC, -1=single
    fImpl->evtTree->Branch("vtx_x",        &fImpl->e_vtx_x);      // mm
    fImpl->evtTree->Branch("vtx_y",        &fImpl->e_vtx_y);
    fImpl->evtTree->Branch("vtx_z",        &fImpl->e_vtx_z);
    fImpl->evtTree->Branch("inv_mass",     &fImpl->e_inv_mass);    // MeV
    fImpl->evtTree->Branch("em_ke",        &fImpl->e_em_ke);
    fImpl->evtTree->Branch("ep_ke",        &fImpl->e_ep_ke);
    fImpl->evtTree->Branch("em_px",        &fImpl->e_em_px);
    fImpl->evtTree->Branch("em_py",        &fImpl->e_em_py);
    fImpl->evtTree->Branch("em_pz",        &fImpl->e_em_pz);
    fImpl->evtTree->Branch("ep_px",        &fImpl->e_ep_px);
    fImpl->evtTree->Branch("ep_py",        &fImpl->e_ep_py);
    fImpl->evtTree->Branch("ep_pz",        &fImpl->e_ep_pz);
    fImpl->evtTree->Branch("openingAngle", &fImpl->e_openingAngle);
    fImpl->evtTree->Branch("neutron_E_eV", &fImpl->e_neutron_E_eV);
    fImpl->evtTree->Branch("capture_vol",  fImpl->e_capture_vol, "capture_vol[32]/C");
    fImpl->evtTree->Branch("cap_x",        &fImpl->e_cap_x);     // mm
    fImpl->evtTree->Branch("cap_y",        &fImpl->e_cap_y);
    fImpl->evtTree->Branch("cap_z",        &fImpl->e_cap_z);

    G4cout << "RunAction: Opened " << fname << G4endl;

#else
    fImpl->hitFile.open(ss.str() + "_hits.csv");
    fImpl->hitFile << "eventID,trackID,parentID,armID,detType,particle,"
                      "u_mm,v_mm,w_mm,edep_eV,ke_MeV,time_ns,"
                      "gx_mm,gy_mm,gz_mm,px,py,pz\n";

    fImpl->evtFile.open(ss.str() + "_events.csv");
    fImpl->evtFile << "eventID,event_type,vtx_x_mm,vtx_y_mm,vtx_z_mm,inv_mass_MeV,"
                      "em_ke_MeV,ep_ke_MeV,"
                      "em_px,em_py,em_pz,ep_px,ep_py,ep_pz,openingAngle_deg\n";

    G4cout << "RunAction: Writing " << ss.str() << "_hits.csv + _events.csv" << G4endl;
#endif
}

// ─────────────────────────────────────────────────────────────────────────────
void RunAction::EndOfRunAction(const G4Run* run) {
    if (!fIsMaster) {
#ifdef USE_ROOT
        if (fImpl->rootFile) {
            fImpl->rootFile->Write();
            fImpl->rootFile->Close();
            delete fImpl->rootFile;
            fImpl->rootFile = nullptr;
        }
#else
        fImpl->hitFile.close();
        fImpl->evtFile.close();
#endif
    }

    // Worker threads print their own per-thread summary.
    // Master prints the G4Run total (which aggregates across all workers in MT).
    if (!fIsMaster) {
#ifdef USE_ROOT
        if (fTotalEvents > 0)
            G4cout << "G4WT" << G4Threading::G4GetThreadId()
                   << " > Wrote " << fTotalEvents << " events, "
                   << fTotalHits  << " hits"
                   << "  (<hits/event>=" << G4double(fTotalHits)/G4double(fTotalEvents) << ")\n";
#endif
    }

    if (fIsMaster || !G4Threading::IsMultithreadedApplication()) {
        G4int nev = run->GetNumberOfEvent();
        G4cout << "\n========= Run Summary =========\n";
        G4cout << "  Events : " << nev << "\n";
        G4cout << "================================\n";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
void RunAction::RecordEvent(const EventData& data) {
    ++fTotalEvents;
    fTotalHits += static_cast<long>(data.hits.size());

#ifdef USE_ROOT
    if (!fImpl->evtTree) return;

    // EventTree
    fImpl->e_eventID      = data.eventID;
    fImpl->e_event_type   = data.event_type;
    fImpl->e_vtx_x        = data.kin.vx;
    fImpl->e_vtx_y        = data.kin.vy;
    fImpl->e_vtx_z        = data.kin.vz;
    fImpl->e_inv_mass     = data.kin.inv_mass_MeV;
    fImpl->e_em_ke        = data.kin.em_ke;
    fImpl->e_ep_ke        = data.kin.ep_ke;
    fImpl->e_em_px        = data.kin.em_px; fImpl->e_em_py = data.kin.em_py; fImpl->e_em_pz = data.kin.em_pz;
    fImpl->e_ep_px        = data.kin.ep_px; fImpl->e_ep_py = data.kin.ep_py; fImpl->e_ep_pz = data.kin.ep_pz;
    fImpl->e_openingAngle = data.kin.openingAngle_deg;
    fImpl->e_neutron_E_eV = data.neutron_E_eV;
    std::strncpy(fImpl->e_capture_vol, data.capture_vol.c_str(), 31);
    fImpl->e_capture_vol[31] = '\0';
    fImpl->e_cap_x = data.cap_x;
    fImpl->e_cap_y = data.cap_y;
    fImpl->e_cap_z = data.cap_z;
    fImpl->evtTree->Fill();

    // HitTree
    for (const auto& h : data.hits) {
        fImpl->h_eventID  = h.eventID;
        fImpl->h_trackID  = h.trackID;
        fImpl->h_parentID = h.parentID;
        fImpl->h_armID    = h.armID;
        std::strncpy(fImpl->h_detType,  h.detType,  31); fImpl->h_detType[31]  = '\0';
        std::strncpy(fImpl->h_particle, h.particle, 31); fImpl->h_particle[31] = '\0';
        fImpl->h_u    = h.u;   fImpl->h_v = h.v;   fImpl->h_w = h.w;
        fImpl->h_edep = h.edep;
        fImpl->h_ke   = h.ke;
        fImpl->h_time = h.time;
        fImpl->h_gx   = h.gx; fImpl->h_gy = h.gy; fImpl->h_gz = h.gz;
        fImpl->h_px   = h.px; fImpl->h_py = h.py; fImpl->h_pz = h.pz;
        fImpl->hitTree->Fill();
    }

#else
    fImpl->evtFile << data.eventID << "," << data.event_type << ","
                   << data.kin.vx << "," << data.kin.vy << "," << data.kin.vz << ","
                   << data.kin.inv_mass_MeV << ","
                   << data.kin.em_ke << "," << data.kin.ep_ke << ","
                   << data.kin.em_px << "," << data.kin.em_py << "," << data.kin.em_pz << ","
                   << data.kin.ep_px << "," << data.kin.ep_py << "," << data.kin.ep_pz << ","
                   << data.kin.openingAngle_deg << "\n";

    for (const auto& h : data.hits) {
        fImpl->hitFile << h.eventID  << "," << h.trackID  << "," << h.parentID << ","
                       << h.armID    << "," << h.detType  << "," << h.particle  << ","
                       << h.u        << "," << h.v        << "," << h.w         << ","
                       << h.edep     << "," << h.ke       << "," << h.time      << ","
                       << h.gx       << "," << h.gy       << "," << h.gz        << ","
                       << h.px       << "," << h.py       << "," << h.pz        << "\n";
    }
#endif
}
