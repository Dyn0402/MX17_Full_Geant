#pragma once
// HitData.hh — one entry in the flat per-hit ROOT tree

#include <cstring>

struct HitData {
    int    eventID  = -1;
    int    trackID  = -1;
    int    parentID = -1;
    int    armID    = -1;   // 0=+X  1=-X  2=+Z  3=-Z
    char   detType[32]  = {};  // "DriftGas", "PlasticScint", "LiqScint_1", …
    char   particle[32] = {};  // "e-", "e+", "gamma", …
    double u    = 0.0;  // local: horizontal-transverse [mm]
    double v    = 0.0;  // local: along beam (Y)       [mm]
    double w    = 0.0;  // local: depth into volume     [mm]
    double edep = 0.0;  // energy deposited             [eV]
    double ke   = 0.0;  // kinetic energy at step start [MeV]
    double time = 0.0;  // global time                  [ns]
    double gx   = 0.0;  // global hit position          [mm]
    double gy   = 0.0;
    double gz   = 0.0;
    double px   = 0.0;  // momentum unit vector
    double py   = 0.0;
    double pz   = 0.0;

    void Clear() { std::memset(this, 0, sizeof(*this)); eventID = armID = trackID = parentID = -1; }
};
