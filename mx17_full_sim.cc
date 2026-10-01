// mx17_full_sim.cc
// MX17 full 4-arm Geant4 simulation.
// X17 → e+e- pairs tracked through the full detector stack.

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "SimConfig.hh"

#include <iostream>
#include <string>
#include <ctime>
#include <algorithm>

static void PrintUsage() {
    std::cerr << "Usage: mx17_full_sim [options] [macro_file]\n"
              << "Options:\n"
              << "  -n <nevents>     Number of events (default: 10000)\n"
              << "  -o <output>      Output file base name (default: x17_output)\n"
              << "  -s <seed>        Random seed (default: time-based)\n"
              << "  -t <nthreads>    MT threads (default: 1)\n"
              << "  -g <gas>         Gas: ArCF4, ArIso, HeEth, ArCO2, etc. (default: ArCF4)\n"
              << "  -v               Verbose output\n"
              << "  --single <name> <E_MeV> <theta_deg> <phi_deg>\n"
              << "                   Single-particle mode (e.g. --single e- 8 90 0)\n"
              << "  --neutron <flux.root> <lambda2d.root>\n"
              << "                   Neutron-beam mode: EAR2 flux + radial profile\n"
              << "                   (data/fluxEAR2-Ph3_in_different_units.root,\n"
              << "                    data/lamda2DvsEn_EAR2.root)\n"
              << "  --bias-ncapture <factor>\n"
              << "                   Scale the ³He(n,γ) cross-section in the gas by <factor>\n"
              << "                   (variance reduction for the rare radiative channel; each\n"
              << "                   biased capture carries weight 1/factor). Neutron mode only.\n"
              << "  --no-al          Replace the He3 capsule Al vessel with vacuum (cross-check:\n"
              << "                   does removing Al kill the capture-gamma background).\n"
              << "  --bias-wall <f> / --bias-thick <f> / --bias-air <f>\n"
              << "                   Same nCapture biasing in the cell's thin walls (window, skin,\n"
              << "                   rods) / thick parts (Al ring + cap, 6LiF) / the World air\n"
              << "  --gamma-cut-um <um>\n"
              << "                   Override the gamma production cut (default 100 um).\n"
              << "  --emin <eV>      Neutron sampling window minimum (default: 1e-3)\n"
              << "  --emax <eV>      Neutron sampling window maximum (default: 1000)\n"
              << "  --gamma-source <capture_lib.csv>\n"
              << "                   Biased wall-background mode: capture-cascade gammas\n"
              << "                   from a capture-vertex library (make_capture_library.py)\n"
              << "  --trajdump [N]   Dump per-step neutron(+secondary) trajectories for the\n"
              << "                   first N events (default 20) to <out>_traj[_t<tid>].csv\n"
              << "                   (event displays; use -t 1 + narrow --emin/--emax)\n"
              << "  --mass <MeV>     X17 mass (default: 16.8)\n"
              << "  --energy <MeV>   4He* transition energy for both X17 and IPC (default: 20.58)\n"
              << "  --ipc <frac>     IPC fraction 0..1 (0=all X17, 1=all IPC, default: 0.5)\n"
              << "  --ipc-multipole ansatz|M1|E0|E1\n"
              << "                   IPC kinematics: legacy 1/M isotropic ansatz (default) or\n"
              << "                   the exact Born distribution of one multipole (ipc_born.py)\n"
              << "  --pair-vertex-lib <lib.csv>\n"
              << "                   Sample X17/IPC vertices from a He3Gas capture-position\n"
              << "                   library (make_capture_library.py --gas-lib) instead of\n"
              << "                   uniformly in the gas (thermal self-shielding profile)\n"
              << "  --pair-vertex-vol <name>  library volume to use (default He3Gas)\n"
              << "  --dist <cm>      Arm distance from target (default: 22.0)\n"
              << "\n  ILL (HANDOFF_SIM.md in x17_facility_search/ill):\n"
              << "  --beam ill <spectrum.csv>\n"
              << "                   Neutron mode with the PF1B/H113 reactor beam (λ from the\n"
              << "                   CSV, divergence κλ, guide-exit acceptance by back-projection)\n"
              << "  --beam-radius <mm>   defining aperture radius (default 10)\n"
              << "  --gun-dist <mm>      aperture → cell entrance window (default 300)\n"
              << "  --exit-dist <mm>     guide exit → aperture (default 1200)\n"
              << "  --vertical-axis x|z  sim axis carrying the exit's 200 mm side (default z)\n"
              << "  --lambda <A>         mono-wavelength instead of the spectrum\n"
              << "  --kappa <rad/A>      divergence slope (default 0.0017; 0 = pencil)\n"
              << "  --slab <Mat:mm>      bare Ø100 mm slab at the origin instead of a target\n"
              << "  --ts | --no-ts       force thermal scattering for solids on/off\n"
              << "                       (default: on with --beam ill, --target cell, --slab)\n"
              << "  --cosmic             cosmic muons from a 3 x 3 m plane 1.5 m above the target\n"
              << "                       (zenith along --vertical-axis; event_type 4)\n"
              << "  --target capsule|cell\n"
              << "  --cell-pressure <bar> --cell-length <mm> --cell-radius <mm> --cell-yw <mm>\n"
              << "  --skin <Mat:mm>  --rods <N>  --window <Mat:mm>  --aperture <mm>\n"
              << "  --end-cap <Mat:mm[+Mat:mm]>  --scraper <rin_mm:t_mm | none>\n"
              << "  --endcap-ring <t_mm:land_mm>\n"
              << "  -h               Print this help\n";
}

int main(int argc, char** argv) {
    SimConfig config;
    config.seed = static_cast<long>(std::time(nullptr));
    std::string macroFile;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if      (a == "-h")              { PrintUsage(); return 0; }
        else if (a == "-n" && i+1<argc)  config.nEvents   = std::stoi(argv[++i]);
        else if (a == "-o" && i+1<argc)  config.outFile   = argv[++i];
        else if (a == "-s" && i+1<argc)  config.seed      = std::stol(argv[++i]);
        else if (a == "-t" && i+1<argc)  config.nThreads  = std::stoi(argv[++i]);
        else if (a == "-g" && i+1<argc)  config.gas       = argv[++i];
        else if (a == "-v")              config.verbose   = true;
        else if (a == "--mass"   && i+1<argc) config.x17Mass_MeV           = std::stod(argv[++i]);
        else if (a == "--energy" && i+1<argc) config.transition_energy_MeV  = std::stod(argv[++i]);
        else if (a == "--ipc"    && i+1<argc) config.ipc_fraction            = std::stod(argv[++i]);
        else if (a == "--pair-vertex-lib" && i+1<argc) config.pairVertexLibFile = argv[++i];
        else if (a == "--pair-vertex-vol" && i+1<argc) config.pairVertexVol = argv[++i];
        else if (a == "--ipc-multipole" && i+1<argc) {
            config.ipcMultipole = argv[++i];
            if (config.ipcMultipole != "ansatz" && config.ipcMultipole != "M1" &&
                config.ipcMultipole != "E0" && config.ipcMultipole != "E1") {
                std::cerr << "--ipc-multipole must be ansatz, M1, E0 or E1\n"; return 1;
            }
        }
        else if (a == "--dist"   && i+1<argc) {   // uniform override of both MM front-face distances
            config.mm_distance_x_cm = config.mm_distance_z_cm = std::stod(argv[++i]);
        }
        else if (a == "--single" && i+4<argc) {
            config.singleParticle           = true;
            config.singleParticleName       = argv[++i];
            config.singleParticleEnergy_MeV = std::stod(argv[++i]);
            config.singleParticleTheta_deg  = std::stod(argv[++i]);
            config.singleParticlePhi_deg    = std::stod(argv[++i]);
        }
        else if (a == "--neutron" && i+2<argc) {
            config.neutronMode        = true;
            config.neutronFluxFile    = argv[++i];
            config.neutronProfileFile = argv[++i];
        }
        else if (a == "--bias-ncapture" && i+1<argc) config.biasNCaptureFactor = std::stod(argv[++i]);
        else if (a == "--bias-wall" && i+1<argc) config.biasWallFactor = std::stod(argv[++i]);
        else if (a == "--bias-air"  && i+1<argc) config.biasAirFactor  = std::stod(argv[++i]);
        else if (a == "--bias-thick" && i+1<argc) config.biasThickFactor = std::stod(argv[++i]);
        else if (a == "--no-al")                     config.disableAlCapsule  = true;
        else if (a == "--gamma-cut-um" && i+1<argc)  config.gammaCut_um       = std::stod(argv[++i]);
        else if (a == "--emin" && i+1<argc) config.neutronEmin_eV = std::stod(argv[++i]);
        else if (a == "--emax" && i+1<argc) config.neutronEmax_eV = std::stod(argv[++i]);
        else if (a == "--gamma-source" && i+1<argc) {
            config.gammaSourceMode = true;
            config.captureLibFile  = argv[++i];
        }
        else if (a == "--trajdump") {
            config.trajDump = true;
            if (i+1 < argc && argv[i+1][0] != '-')
                config.trajDumpMaxEvents = std::stoi(argv[++i]);
        }
        else if (a == "--beam" && i+2<argc) {
            std::string kind = argv[++i];
            if (kind != "ill") { std::cerr << "Unknown beam: " << kind << "\n"; return 1; }
            config.neutronMode     = true;
            config.illBeam         = true;
            config.illSpectrumFile = argv[++i];
        }
        else if (a == "--beam-radius" && i+1<argc) config.illBeamRadius_mm = std::stod(argv[++i]);
        else if (a == "--gun-dist"    && i+1<argc) config.illGunDist_mm    = std::stod(argv[++i]);
        else if (a == "--exit-dist"   && i+1<argc) config.illExitDist_mm   = std::stod(argv[++i]);
        else if (a == "--vertical-axis" && i+1<argc) {
            std::string ax = argv[++i];
            if (ax != "x" && ax != "z") { std::cerr << "--vertical-axis must be x or z\n"; return 1; }
            config.illVerticalAxis = ax[0];
        }
        else if (a == "--lambda" && i+1<argc) config.illLambdaFixed_A   = std::stod(argv[++i]);
        else if (a == "--kappa"  && i+1<argc) config.illKappa_rad_per_A = std::stod(argv[++i]);
        else if (a == "--slab"   && i+1<argc) config.slab               = argv[++i];
        else if (a == "--cosmic") config.cosmic = true;
        else if (a == "--ts")    config.thermalScattering = 1;
        else if (a == "--no-ts") config.thermalScattering = 0;
        else if (a == "--target" && i+1<argc) {
            std::string t = argv[++i];
            if      (t == "cell")    config.cellTarget = true;
            else if (t == "capsule") config.cellTarget = false;
            else { std::cerr << "Unknown target: " << t << "\n"; return 1; }
        }
        else if (a == "--cell-pressure" && i+1<argc) config.cellPressure_bar = std::stod(argv[++i]);
        else if (a == "--cell-length"   && i+1<argc) config.cellLength_mm    = std::stod(argv[++i]);
        else if (a == "--cell-radius"   && i+1<argc) config.cellRadius_mm    = std::stod(argv[++i]);
        else if (a == "--cell-yw"       && i+1<argc) config.cellYw_mm        = std::stod(argv[++i]);
        else if (a == "--skin"          && i+1<argc) config.cellSkin         = argv[++i];
        else if (a == "--rods"          && i+1<argc) config.cellRods         = std::stoi(argv[++i]);
        else if (a == "--window"        && i+1<argc) config.cellWindow       = argv[++i];
        else if (a == "--aperture"      && i+1<argc) config.cellApertureR_mm = std::stod(argv[++i]);
        else if (a == "--end-cap"       && i+1<argc) config.cellEndCap       = argv[++i];
        else if (a == "--scraper"       && i+1<argc) {
            std::string s = argv[++i];
            if (s == "none") config.cellScraperRin_mm = 0.0;
            else {
                auto c = s.find(':');
                config.cellScraperRin_mm = std::stod(s.substr(0, c));
                if (c != std::string::npos) config.cellScraperT_mm = std::stod(s.substr(c + 1));
            }
        }
        else if (a == "--endcap-ring"   && i+1<argc) {
            std::string s = argv[++i];
            auto c = s.find(':');
            config.cellRingThick_mm = std::stod(s.substr(0, c));
            if (c != std::string::npos) config.cellRingLand_mm = std::stod(s.substr(c + 1));
        }
        else if (a[0] != '-') macroFile = a;
        else { std::cerr << "Unknown option: " << a << "\n"; PrintUsage(); return 1; }
    }

    if (config.thermalScattering < 0)
        config.thermalScattering = (config.illBeam || config.cellTarget || !config.slab.empty()) ? 1 : 0;

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(config.seed);

    std::cout << "=== MX17 Full Simulation ===\n"
              << "  Geant4   : " << G4Version << "\n"
              << "  Gas      : " << config.gas << "\n"
              << "  Events   : " << config.nEvents << "\n"
              << "  Output   : " << config.outFile << "\n"
              << "  Seed     : " << config.seed << "\n"
              << "  Threads  : " << config.nThreads << "\n"
              << "  Al vessel: " << (config.disableAlCapsule ? "DISABLED (vacuum)" : "enabled") << "\n"
              << "  Gamma cut: " << config.gammaCut_um << " um\n"
              << "  Therm.sc.: " << (config.thermalScattering ? "on" : "off") << "\n";
    if (config.cellTarget)
        std::cout << "  Target   : 3He cell " << config.cellPressure_bar << " bar, L="
                  << config.cellLength_mm << " mm, R=" << config.cellRadius_mm
                  << " mm, y_w=" << config.cellYw_mm << " mm, skin " << config.cellSkin
                  << ", rods " << config.cellRods << ", window " << config.cellWindow
                  << ", end cap " << config.cellEndCap << "\n";
    if (config.cosmic) {
        const double A_cm2 = config.cosmicPlane_mm * config.cosmicPlane_mm / 100.0;
        std::cout << "  Mode     : cosmic muons, plane " << config.cosmicPlane_mm << " mm at "
                  << config.cosmicHeight_mm << " mm along +" << config.illVerticalAxis
                  << "; live time = N x " << 60.0 / A_cm2 << " s\n";
    }
    else if (config.illBeam)
        std::cout << "  Mode     : ILL PF1B beam  spectrum=" << config.illSpectrumFile
                  << "  aperture r=" << config.illBeamRadius_mm << " mm at "
                  << config.illGunDist_mm << " mm upstream of the window, exit "
                  << config.illExitDist_mm << " mm upstream; vertical="
                  << config.illVerticalAxis << "\n";
    else if (config.neutronMode)
        std::cout << "  Mode     : neutron beam  E=[" << config.neutronEmin_eV
                  << ", " << config.neutronEmax_eV << "] eV\n"
                  << "  Flux     : " << config.neutronFluxFile << "\n"
                  << "  Profile  : " << config.neutronProfileFile << "\n";
    else if (config.gammaSourceMode)
        std::cout << "  Mode     : gamma-source (biased wall background)\n"
                  << "  Library  : " << config.captureLibFile << "\n";
    else if (config.singleParticle)
        std::cout << "  Mode     : single-particle " << config.singleParticleName
                  << " " << config.singleParticleEnergy_MeV << " MeV\n";
    else {
        std::cout << "  Mode     : X17+IPC pairs  m_X17=" << config.x17Mass_MeV
                  << " MeV  E_transition=" << config.transition_energy_MeV
                  << " MeV  ipc_fraction=" << config.ipc_fraction
                  << "  ipc=" << config.ipcMultipole << "\n";
        if (!config.pairVertexLibFile.empty())
            std::cout << "  Vertices : " << config.pairVertexLibFile << "\n";
    }
    std::cout << "============================\n";

#ifdef G4MULTITHREADED
    auto* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(config.nThreads);
#else
    auto* runManager = new G4RunManager;
#endif

    auto* detCon = new DetectorConstruction(config);
    runManager->SetUserInitialization(detCon);
    runManager->SetUserInitialization(new PhysicsList(
        std::max({config.biasNCaptureFactor, config.biasWallFactor, config.biasAirFactor,
                  config.biasThickFactor}),
        config.gammaCut_um,
                                                      config.thermalScattering == 1));
    runManager->SetUserInitialization(new ActionInitialization(config, detCon));
    runManager->Initialize();

    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (macroFile.empty())
        UI->ApplyCommand("/run/beamOn " + std::to_string(config.nEvents));
    else
        UI->ApplyCommand("/control/execute " + macroFile);

    delete runManager;
    if (config.cosmic)
        std::cout << "Cosmic live time: " << config.nEvents * 60.0 /
                         (config.cosmicPlane_mm * config.cosmicPlane_mm / 100.0)
                  << " s\n";
    std::cout << "Done. Output: " << config.outFile << "\n";
    return 0;
}
