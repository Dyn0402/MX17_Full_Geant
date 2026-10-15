#pragma once
// GasProperties.hh
//
// One table for every counting gas the simulation can use. Both the G4Material
// (DetectorConstruction) and the edep -> ion-pair conversion (SteppingAction)
// read from it, so the density, composition, W and Fano factor of a gas cannot
// drift apart again.
//
// Fixed 2026-10-15 (Heed cross-check, MX17_Geant response/validation/heed_vs_g4/;
// ported from MX17_Geant e17c72c):
//  - Densities were 0 C STP values on materials declared at 293.15 K, so every
//    gas was 7.3 % too dense, giving 7.3 % too many primaries.
//  - Mixture fractions are VOLUME fractions but were passed to
//    G4Material::AddMaterial, which takes MASS fractions.
//  - W came from a hand-written per-mixture table (NeIso 27 eV where the
//    components give 35 eV; ArCF4Iso 33 eV where they give 27).
// W and F now follow Garfield++ (MediumGas::GetGasInfo: ICRU 31 W values, IAEA
// TECDOC 799 Fano factors), with the mixture value linear in the mole fraction
// as Garfield does it, so the counts match what Heed is configured with.

#include <string>
#include <utility>
#include <vector>

namespace gas {

// Conditions every gas material is declared at.
constexpr double kTemperatureK = 293.15;
constexpr double kPressurePa   = 101325.0;

struct Component {
    std::string name;                                   // "Ar", "iC4H10", ...
    std::string material;                               // G4Material name
    double molarMass;                                   // g/mol
    std::vector<std::pair<std::string, int>> atoms;     // element symbol, count
    double W;                                           // eV per ion pair
    double fano;
};

struct Mixture {
    std::string name;                                   // SimConfig::gas key = G4Material name
    std::vector<std::pair<std::string, double>> parts;  // component, volume fraction
};

const std::vector<Component>& Components();
const std::vector<Mixture>&   Mixtures();
const Component& FindComponent(const std::string& name);
const Mixture&   FindMixture(const std::string& name);   // throws if unknown

double ComponentDensity(const std::string& component);   // g/cm3, ideal gas
double MixtureDensity(const Mixture& mix);                // g/cm3
double MassFraction(const Mixture& mix, std::size_t i);
double MixtureW(const std::string& name);                 // eV
double MixtureFano(const std::string& name);

}  // namespace gas
