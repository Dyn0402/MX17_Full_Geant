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
              << "  --mass <MeV>     X17 mass (default: 16.8)\n"
              << "  --energy <MeV>   4He* transition energy for both X17 and IPC (default: 20.58)\n"
              << "  --ipc <frac>     IPC fraction 0..1 (0=all X17, 1=all IPC, default: 0.5)\n"
              << "  --dist <cm>      Arm distance from target (default: 22.0)\n"
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
        else if (a == "--dist"   && i+1<argc) config.mm_distance_cm          = std::stod(argv[++i]);
        else if (a == "--single" && i+4<argc) {
            config.singleParticle           = true;
            config.singleParticleName       = argv[++i];
            config.singleParticleEnergy_MeV = std::stod(argv[++i]);
            config.singleParticleTheta_deg  = std::stod(argv[++i]);
            config.singleParticlePhi_deg    = std::stod(argv[++i]);
        }
        else if (a[0] != '-') macroFile = a;
        else { std::cerr << "Unknown option: " << a << "\n"; PrintUsage(); return 1; }
    }

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(config.seed);

    std::cout << "=== MX17 Full Simulation ===\n"
              << "  Geant4   : " << G4Version << "\n"
              << "  Gas      : " << config.gas << "\n"
              << "  Events   : " << config.nEvents << "\n"
              << "  Output   : " << config.outFile << "\n"
              << "  Seed     : " << config.seed << "\n"
              << "  Threads  : " << config.nThreads << "\n";
    if (config.singleParticle)
        std::cout << "  Mode     : single-particle " << config.singleParticleName
                  << " " << config.singleParticleEnergy_MeV << " MeV\n";
    else
        std::cout << "  Mode     : X17+IPC pairs  m_X17=" << config.x17Mass_MeV
                  << " MeV  E_transition=" << config.transition_energy_MeV
                  << " MeV  ipc_fraction=" << config.ipc_fraction << "\n";
    std::cout << "============================\n";

#ifdef G4MULTITHREADED
    auto* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(config.nThreads);
#else
    auto* runManager = new G4RunManager;
#endif

    auto* detCon = new DetectorConstruction(config);
    runManager->SetUserInitialization(detCon);
    runManager->SetUserInitialization(new PhysicsList());
    runManager->SetUserInitialization(new ActionInitialization(config, detCon));
    runManager->Initialize();

    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (macroFile.empty())
        UI->ApplyCommand("/run/beamOn " + std::to_string(config.nEvents));
    else
        UI->ApplyCommand("/control/execute " + macroFile);

    delete runManager;
    std::cout << "Done. Output: " << config.outFile << "\n";
    return 0;
}
