// GasProperties.cc -- see GasProperties.hh for why this table exists.

#include "GasProperties.hh"

#include <cmath>
#include <stdexcept>

namespace gas {

// W [eV] and Fano factor from Garfield++ MediumGas::GetGasInfo (ICRU 31 W values;
// Fano from IAEA TECDOC 799, or the Krajcar-Bronic relation where unmeasured).
const std::vector<Component>& Components() {
    static const std::vector<Component> c = {
        {"Ar",     "PureArgon",  39.948,    {{"Ar", 1}},          26.4, 0.17},
        {"Ne",     "PureNe",     20.1797,   {{"Ne", 1}},          35.4, 0.17},
        {"He",     "PureHe",      4.002602, {{"He", 1}},          41.3, 0.17},
        {"CO2",    "CO2_gas",    44.0095,   {{"C", 1}, {"O", 2}}, 33.0, 0.32},
        {"CF4",    "CF4_gas",    88.0043,   {{"C", 1}, {"F", 4}}, 34.3, 0.26},
        {"C2H6",   "Ethane",     30.069,    {{"C", 2}, {"H", 6}}, 25.0, 0.28},
        {"iC4H10", "Isobutane",  58.1222,   {{"C", 4}, {"H", 10}},23.4, 0.26},
    };
    return c;
}

const std::vector<Mixture>& Mixtures() {
    static const std::vector<Mixture> m = {
        {"ArCF4", {{"Ar", 0.90}, {"CF4", 0.10}}},
        {"HeEth", {{"He", 0.965}, {"C2H6", 0.035}}},
        {"ArCO2", {{"Ar", 0.70}, {"CO2", 0.30}}},
        {"ArIso", {{"Ar", 0.95}, {"iC4H10", 0.05}}},
        {"NeIso", {{"Ne", 0.95}, {"iC4H10", 0.05}}},
        {"NeCF4", {{"Ne", 0.90}, {"CF4", 0.10}}},
        {"ArCF4Iso", {{"Ar", 0.88}, {"CF4", 0.10}, {"iC4H10", 0.02}}},
        {"ArCF4CO2", {{"Ar", 0.45}, {"CF4", 0.40}, {"CO2", 0.15}}},
        {"PureCF4", {{"CF4", 1.0}}},
        {"PureAr", {{"Ar", 1.0}}},
        {"PureHe", {{"He", 1.0}}},
        {"PureNe", {{"Ne", 1.0}}},
        {"PureEthane", {{"C2H6", 1.0}}},
        {"PureIso", {{"iC4H10", 1.0}}},
        {"PureCO2", {{"CO2", 1.0}}},
    };
    return m;
}

const Component& FindComponent(const std::string& name) {
    for (const auto& c : Components())
        if (c.name == name) return c;
    throw std::runtime_error("gas: unknown component " + name);
}

const Mixture& FindMixture(const std::string& name) {
    for (const auto& m : Mixtures())
        if (m.name == name) return m;
    throw std::runtime_error("gas: unknown gas " + name);
}

double ComponentDensity(const std::string& component) {
    // Ideal gas at the declared conditions. Good to <1 % for the noble gases,
    // CO2 and CF4 at 1 atm; isobutane is ~2 % denser than ideal at 20 C, which
    // is <0.1 % on a mixture with 5 % of it.
    const double R = 8.314462618;  // J/(mol K)
    return kPressurePa * FindComponent(component).molarMass * 1e-3 /
           (R * kTemperatureK) * 1e-3;  // kg/m3 -> g/cm3
}

double MixtureDensity(const Mixture& mix) {
    // Partial densities of ideal gases at a common T and P add.
    double fsum = 0.0, rho = 0.0;
    for (const auto& [c, f] : mix.parts) {
        fsum += f;
        rho  += f * ComponentDensity(c);
    }
    if (std::abs(fsum - 1.0) > 1e-6)
        throw std::runtime_error("gas " + mix.name + ": volume fractions must sum to 1");
    return rho;
}

double MassFraction(const Mixture& mix, std::size_t i) {
    const auto& [c, f] = mix.parts.at(i);
    return f * ComponentDensity(c) / MixtureDensity(mix);
}

// Linear in mole fraction, as Garfield (MediumGas::SetComposition) feeds Heed.
double MixtureW(const std::string& name) {
    double w = 0.0;
    for (const auto& [c, f] : FindMixture(name).parts) w += f * FindComponent(c).W;
    return w;
}

double MixtureFano(const std::string& name) {
    double fano = 0.0;
    for (const auto& [c, f] : FindMixture(name).parts) fano += f * FindComponent(c).fano;
    return fano;
}

}  // namespace gas
