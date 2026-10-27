// OgiveProfile.hh — Pierfrancesco's pointed target shell (his CAD sketch of 27 Oct 2026,
// x17_facility_search lnl/mech_models/sketch li target X17.step; TARGET_OPTIONS.md §19,
// SHAPES_GEANT.md L11).  A thin shell of revolution about the beam axis: a "lemon" flank
// (an arc of radius A whose centre sits RM off the axis on the far side) closed by a
// spherical tip.  Exact B-rep numbers in mm, in a frame where the INNER tip sits at y = 0
// (the CAD frame shifted by −2.993 mm) and the beam comes along +y into the open base.
//
//   outer flank arc A_OUT = 375.286576, inner A_IN = 375.086576 (wall 0.200), both centred at
//   (r = −RM, y = YC); tip sphere centre (0, TIP_C, 0), radii TIP_RO = 3.0 / TIP_RI = 2.8;
//   the flank meets the tip at y = Y_CAP_O / Y_CAP_I; open base at y = Y_BASE (Ø50 outer).
//
// The same functions build the G4Polycone layers (DetectorConstruction) and place the
// pair / γ vertices in the film (X17PrimaryGenerator), so the two cannot disagree.
// Layers inside the shell (a liner, the film) are the inner surface offset inward along
// its normal: for a circle arc that is the same centre with a smaller radius, for the
// sphere a smaller radius, and the flank/tip junction moves along the junction normal.
#ifndef OGIVE_PROFILE_HH
#define OGIVE_PROFILE_HH

#include <algorithm>
#include <cmath>
#include <vector>

namespace ogive {

constexpr double Y_SHIFT = -2.993;                    // CAD y → this frame
constexpr double RM      = 337.656185;                // lemon centre, off-axis distance (mm)
constexpr double YC      = -156.606553 + Y_SHIFT;     // lemon centre y
constexpr double A_OUT   = 375.286576;                // flank arc radii (outer, inner)
constexpr double A_IN    = 375.086576;
constexpr double TIP_C   = 0.193096 + Y_SHIFT;        // tip sphere centre y
constexpr double TIP_RO  = 3.0;
constexpr double TIP_RI  = 2.8;
constexpr double Y_CAP_O = 1.5   + Y_SHIFT;           // flank / tip junction (outer, inner)
constexpr double Y_CAP_I = 1.414 + Y_SHIFT;
constexpr double Y_BASE  = -60.0 + Y_SHIFT;           // open base
constexpr double WALL    = 0.200;                     // mm

// One surface of revolution: the inner wall offset inward by `off` mm (off = 0 is the
// shell's inner surface; off = −WALL its outer surface).
struct Surface {
    double a, rt, ycap;                               // flank radius, tip radius, junction y
    explicit Surface(double off_mm = 0.0)
        : a(A_IN - off_mm), rt(TIP_RI - off_mm),
          ycap(TIP_C + (Y_CAP_I - TIP_C) * (TIP_RI - off_mm) / TIP_RI) {}
    double yTop() const { return TIP_C + rt; }
    // radius at height y (0 above the tip, 0 below the base)
    double r(double y) const {
        if (y < Y_BASE || y > yTop()) return 0.0;
        double rr = 0.0;
        if (y <= ycap) {
            const double s = a * a - (y - YC) * (y - YC);
            if (s > 0) rr = std::max(rr, std::sqrt(s) - RM);
        }
        const double t = rt * rt - (y - TIP_C) * (y - TIP_C);
        if (t > 0) rr = std::max(rr, std::sqrt(t));
        return std::max(rr, 0.0);
    }
    // the surface point above (x, z) = (rho) seen along +y: its y and the inward unit normal
    // (nr radial component, ny axial).  Returns false if rho is outside the base.
    bool hit(double rho, double& y, double& nr, double& ny) const {
        double yc = -1e9, yf = -1e9;
        const double t = rt * rt - rho * rho;
        if (t >= 0) yc = TIP_C + std::sqrt(t);
        const double s = a * a - (rho + RM) * (rho + RM);
        if (s >= 0) { yf = YC + std::sqrt(s); if (yf > ycap) yf = -1e9; }
        if (yc < -1e8 && yf < -1e8) return false;
        if (yc >= yf) {                                // on the tip sphere
            y = yc;
            const double n = std::hypot(rho, y - TIP_C);
            nr = -rho / n; ny = -(y - TIP_C) / n;
        } else {                                       // on the flank
            y = yf;
            const double n = std::hypot(rho + RM, y - YC);
            nr = -(rho + RM) / n; ny = -(y - YC) / n;
        }
        return y >= Y_BASE;
    }
};

// z-planes for the polycone layers: 0.5 mm along the flank, 0.05 mm over the tip region
inline std::vector<double> planes(double yTop) {
    std::vector<double> z;
    const double yFine = TIP_C - TIP_RO - 0.5;
    for (double y = Y_BASE; y < yFine; y += 0.5) z.push_back(y);
    for (double y = yFine; y < yTop - 1e-6; y += 0.05) z.push_back(y);
    z.push_back(yTop);
    return z;
}

}   // namespace ogive
#endif
