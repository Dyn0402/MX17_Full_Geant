// DetectorConstruction.cc
// 4-arm X17 experiment geometry.
// Beam axis = +Y. Arms at ±X (arm 0,1) and ±Z (arm 2,3).
//
// Stack per arm, inside → out (measured 2026-07-15; depths from the MM
// drift-mylar front face; see GEOMETRY_CHANGE_CHECKLIST.md):
//   MM window layers | MM drift gas (30 mm, scored) | amp gas | resistive paste
//   PCB stack
//   [air]  — gap grows so the SiPM container front lands 11 cm from mylar front
//   SiPM trigger wall: 16 of 20 plastic scint bars (2.5 cm each), scored as
//     "PlasticScint", centered on the STRUCTURE (u=0), read-out window shifted
//     1 bar toward the MM.  3 mm scint centered in a 3.3 cm container.
//   [air, per-arm gap measured 2026-07-17]
//   Plastics: 2× (black-mylar tape envelope | 20×30×2.5cm plastic scint (scored))
//     side-by-side in u, centered on the MM (pinwheel-shifted).
//   [air]
//   LS vessel (STEP "LS X17.step"): 2.6 mm CFRP shell — 45×45 cm slab +
//     funnel + Ø50 mm neck with half-inserted PMT; LAB liquid (scored) fills
//     slab+funnel+neck and bulges the slab faces (6.5 L fill).  Surveyed
//     2026-07-17: per-arm depth (flat slab face from SiPM back), per-arm
//     height offset, VERTICAL (PMT up) on B/C, HORIZONTAL (PMT +u) on A/D.

#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Ellipsoid.hh"
#include "G4Polycone.hh"
#include "G4UnionSolid.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"

#include <stdexcept>
#include <cmath>

DetectorConstruction::DetectorConstruction(const SimConfig& cfg)
    : G4VUserDetectorConstruction(), fConfig(cfg) {}

// ─────────────────────────────────────────────────────────────────────────────
void DetectorConstruction::DefineMaterials() {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* elH  = nist->FindOrBuildElement("H");
    G4Element* elC  = nist->FindOrBuildElement("C");
    G4Element* elN  = nist->FindOrBuildElement("N");
    G4Element* elO  = nist->FindOrBuildElement("O");
    G4Element* elSi = nist->FindOrBuildElement("Si");
    G4Element* elCl = nist->FindOrBuildElement("Cl");
    G4Element* elAr = nist->FindOrBuildElement("Ar");
    G4Element* elNe = nist->FindOrBuildElement("Ne");
    G4Element* elHe = nist->FindOrBuildElement("He");
    G4Element* elF  = nist->FindOrBuildElement("F");

    // ── Gas component materials ──────────────────────────────
    auto* isobutane = new G4Material("Isobutane", 2.67e-3*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
    isobutane->AddElement(elC, 4); isobutane->AddElement(elH, 10);

    auto* ethane = new G4Material("Ethane", 1.356e-3*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
    ethane->AddElement(elC, 2); ethane->AddElement(elH, 6);

    auto* CO2 = new G4Material("CO2_gas", 1.977e-3*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
    CO2->AddElement(elC, 1); CO2->AddElement(elO, 2);

    auto* CF4 = new G4Material("CF4_gas", 3.72e-3*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
    CF4->AddElement(elC, 1); CF4->AddElement(elF, 4);

    auto* purAr = new G4Material("PureArgon", 1.782e-3*g/cm3, 1, kStateGas, 293.15*kelvin, 1*atmosphere);
    purAr->AddElement(elAr, 1);

    auto* pureHe = new G4Material("PureHe", 0.1786e-3*g/cm3, 1, kStateGas, 293.15*kelvin, 1*atmosphere);
    pureHe->AddElement(elHe, 1);

    auto* pureNe = new G4Material("PureNe", 0.8999e-3*g/cm3, 1, kStateGas, 293.15*kelvin, 1*atmosphere);
    pureNe->AddElement(elNe, 1);

    auto makeMix2 = [&](const char* name, G4double rho,
                        G4Material* m1, G4double f1,
                        G4Material* m2, G4double f2) -> G4Material* {
        auto* m = new G4Material(name, rho*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(m1, f1); m->AddMaterial(m2, f2); return m;
    };
    auto makeMix3 = [&](const char* name, G4double rho,
                        G4Material* m1, G4double f1,
                        G4Material* m2, G4double f2,
                        G4Material* m3, G4double f3) -> G4Material* {
        auto* m = new G4Material(name, rho*g/cm3, 3, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(m1, f1); m->AddMaterial(m2, f2); m->AddMaterial(m3, f3); return m;
    };

    fMats["ArCF4"]    = makeMix2("ArCF4",    0.90*1.782e-3+0.10*3.72e-3, purAr,0.90,CF4,0.10);
    fMats["HeEth"]    = makeMix2("HeEth",    0.965*0.1786e-3+0.035*1.356e-3, pureHe,0.965,ethane,0.035);
    fMats["ArCO2"]    = makeMix2("ArCO2",    0.70*1.782e-3+0.30*1.977e-3, purAr,0.70,CO2,0.30);
    fMats["ArIso"]    = makeMix2("ArIso",    0.95*1.782e-3+0.05*2.67e-3, purAr,0.95,isobutane,0.05);
    fMats["NeIso"]    = makeMix2("NeIso",    0.95*0.8999e-3+0.05*2.67e-3, pureNe,0.95,isobutane,0.05);
    fMats["NeCF4"]    = makeMix2("NeCF4",    0.90*0.8999e-3+0.10*3.72e-3, pureNe,0.90,CF4,0.10);
    fMats["ArCF4Iso"] = makeMix3("ArCF4Iso", 0.88*1.782e-3+0.10*3.72e-3+0.02*2.67e-3, purAr,0.88,CF4,0.10,isobutane,0.02);
    fMats["ArCF4CO2"] = makeMix3("ArCF4CO2", 0.45*1.782e-3+0.40*3.72e-3+0.15*1.977e-3, purAr,0.45,CF4,0.40,CO2,0.15);
    auto* purCF4 = new G4Material("PureCF4", 3.72e-3*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
    purCF4->AddElement(elC,1); purCF4->AddElement(elF,4);
    fMats["PureCF4"]   = purCF4;
    fMats["PureAr"]    = purAr;
    fMats["PureHe"]    = pureHe;
    fMats["PureNe"]    = pureNe;
    fMats["PureEthane"]= ethane;
    fMats["PureIso"]   = isobutane;
    fMats["PureCO2"]   = CO2;

    // ── He-3 gas at 500 bar ──────────────────────────────────
    {
        auto* isoHe3 = new G4Isotope("He3_iso", 2, 3, 3.0160293*g/mole);
        auto* elHe3  = new G4Element("Helium3_elem", "3He", 1);
        elHe3->AddIsotope(isoHe3, 1.0);
        auto* m = new G4Material("He3Gas_500bar", 62.7e-3*g/cm3, 1,
                                  kStateGas, 293.15*kelvin, 500*atmosphere);
        m->AddElement(elHe3, 1);
        fMats["He3Gas"] = m;
    }

    // ── Structural / detector materials ─────────────────────
    {
        auto* m = new G4Material("CFRP", 1.55*g/cm3, 3);
        m->AddElement(elC, 0.8968); m->AddElement(elH, 0.0207); m->AddElement(elO, 0.0826);
        fMats["CFRP"] = m;
    }
    {
        auto* m = new G4Material("ResistivePaste", 1.4*g/cm3, 3);
        m->AddElement(elC, 0.65); m->AddElement(elH, 0.08); m->AddElement(elO, 0.27);
        fMats["ResistivePaste"] = m;
    }
    {
        auto* m = new G4Material("FR4", 1.85*g/cm3, 4);
        m->AddElement(elSi, 0.2805); m->AddElement(elO, 0.4195);
        m->AddElement(elC,  0.2750); m->AddElement(elH, 0.0250);
        fMats["FR4"] = m;
    }
    {
        auto* m = new G4Material("Rohacell51", 0.052*g/cm3, 4);
        m->AddElement(elC, 0.5783); m->AddElement(elH, 0.0602);
        m->AddElement(elN, 0.1687); m->AddElement(elO, 0.1928);
        fMats["Rohacell51"] = m;
    }
    {
        auto* m = new G4Material("LAB_LiqScint", 0.86*g/cm3, 2);
        m->AddElement(elC, 0.8780); m->AddElement(elH, 0.1220);
        fMats["LAB"] = m;
    }
    // Black mylar tape: used for trigger scint wall wrapping and back scint wrapping.
    // PET (polyethylene terephthalate), density 1.38 g/cm3.
    // Black dye negligible for radiation; treated as pure PET.
    fMats["BlackMylar"] = nist->FindOrBuildMaterial("G4_MYLAR");
}

G4Material* DetectorConstruction::GetMat(const std::string& name) const {
    auto it = fMats.find(name);
    if (it == fMats.end())
        throw std::runtime_error("DetectorConstruction: unknown material '" + name + "'");
    return it->second;
}

// ─────────────────────────────────────────────────────────────────────────────
G4VPhysicalVolume* DetectorConstruction::Construct() {
    DefineMaterials();

    G4NistManager* nist = G4NistManager::Instance();
    G4Material* matAir    = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* matMylar  = nist->FindOrBuildMaterial("G4_MYLAR");
    G4Material* matAl     = nist->FindOrBuildMaterial("G4_Al");
    G4Material* matKapton = nist->FindOrBuildMaterial("G4_KAPTON");
    G4Material* matCu     = nist->FindOrBuildMaterial("G4_Cu");
    G4Material* matSteel  = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* matGas    = GetMat(fConfig.gas);
    G4Material* matPlScint= nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4Material* matBlkMylar = GetMat("BlackMylar");

    fMats["_Gas"] = matGas;

    // ── Layer thicknesses ────────────────────────────────────
    // MM stack
    G4double tMylar    = 40.0  * um;
    G4double tAlWin    = 0.1   * um;
    G4double tKapCath  = 50.0  * um;
    G4double tCuCath   = 9.0   * um;
    G4double tDrift    = 30.0  * mm;
    G4double tMesh     = 30.0  * um;
    G4double tAmp      = 150.0 * um;
    G4double tResPaste = 100.0 * um;

    // PCB stack
    G4double tPCBKap  = 50.0  * um;
    G4double tPCBCu   = 26.0  * um;
    G4double tPCBFR4  = 100.0 * um;
    G4double tPCBRoh  = 5.0   * mm;
    G4double tPCBAl   = 50.0  * um;

    // SiPM trigger wall: bare plastic scint bar (2.5 cm wide, 50 cm long)
    G4double tSipmScint = fConfig.sipm_scint_thick_cm * cm;  // 3 mm active depth

    // LS vessel (STEP "LS X17.step"; see SimConfig.hh note).  Vessel local
    // frame: x = u (slab width), y = v (long axis, neck at +v), z = w (depth).
    G4double lsUo   = fConfig.ls_slab_u_cm     * cm / 2;  // slab outer half-width
    G4double lsVo   = fConfig.ls_slab_v_cm     * cm / 2;  // slab outer half-length
    G4double lsTo   = fConfig.ls_slab_thick_cm * cm / 2;  // slab outer half-thickness
    G4double lsWall = fConfig.ls_wall_mm       * mm;      // CFRP wall
    G4double lsFunL = fConfig.ls_funnel_len_cm * cm;      // funnel length
    G4double lsNkL  = fConfig.ls_neck_len_cm   * cm;      // neck length
    G4double lsNkR  = fConfig.ls_neck_r_cm     * cm;      // neck outer radius
    // Interior (liquid) dims: outer shrunk by the wall
    G4double lsUi   = lsUo - lsWall;
    G4double lsTi   = lsTo - lsWall;
    G4double lsNkRi = lsNkR - lsWall;
    G4double lsSlabInHV  = lsVo - lsWall/2;   // interior slab half-length (v)
    G4double lsSlabInCen = lsWall/2;          // its centre offset toward the neck
    // PMT: window inserted ls_pmt_insert_frac into the neck
    G4double pmtR     = fConfig.ls_pmt_r_cm     * cm;
    G4double pmtGlass = fConfig.ls_pmt_glass_mm * mm;
    G4double pmtIns   = fConfig.ls_pmt_insert_frac * lsNkL;
    G4double pmtOut   = std::max(0.0, fConfig.ls_pmt_len_cm * cm - pmtIns);
    G4double pmtFaceV = lsVo + lsFunL + pmtIns;   // PMT window plane (from slab centre)
    // Interior funnel: the outer loft shrunk by the wall, truncated one wall
    // thickness before the neck so the funnel's square end keeps CFRP at the
    // corners; a round stub bridges from there to the PMT window.
    G4double stubStartV = lsVo + lsFunL - lsWall;
    G4double funSlopeU  = (lsUi - lsNkRi) / lsFunL;
    G4double funSlopeT  = (lsNkRi - lsTi) / lsFunL;
    G4double iFunHL = (lsFunL - lsWall) / 2;           // interior funnel half-length
    G4double iFunU2 = lsNkRi + funSlopeU * lsWall;     // half-width at its narrow end
    G4double iFunT2 = lsNkRi - funSlopeT * lsWall;     // half-thickness at its narrow end
    // Bulge: the fill volume exceeds the flat interior; the excess is carried
    // by two ellipsoidal domes on the slab faces (height solved from volume).
    G4double vSlabIn = (2*lsUi) * (2*lsSlabInHV) * (2*lsTi);
    G4double fA1 = (2*lsUi)*(2*lsTi), fA2 = (2*iFunU2)*(2*iFunT2);
    G4double fAm = (lsUi + iFunU2) * (lsTi + iFunT2);
    G4double vFunIn  = (2*iFunHL)/6.0 * (fA1 + 4*fAm + fA2);   // prismatoid
    G4double vStubIn = pi * lsNkRi*lsNkRi * (pmtFaceV - stubStartV);
    G4double vFill   = fConfig.ls_fill_liters * 1000.0 * cm3;
    G4double hCap    = std::max(0.5*mm,
        (vFill - vSlabIn - vFunIn - vStubIn) / 2.0 / (2.0/3.0 * pi * lsUi * lsSlabInHV));

    // Plastic scint wrapping: tape (outer) → Al foil → PVT
    G4double tTape    = fConfig.backscint_tape_um * um;  // 200 µm black mylar
    G4double tBscAl   = fConfig.backscint_al_um   * um;  // 20 µm Al foil

    // ── Half-sizes in arm local frame (u, v) ────────────────
    G4double mmU_hf  = fConfig.mm_size_u_cm    * cm / 2;
    G4double mmV_hf  = fConfig.mm_size_v_cm    * cm / 2;

    // SiPM bar half-sizes
    G4double sipmBar_hu = fConfig.sipm_bar_width_cm * cm / 2;  // 1.25 cm
    G4double sipmBar_hv = fConfig.sipm_size_v_cm    * cm / 2;  // 25 cm
    G4double sipmBar_hw = tSipmScint / 2;

    // Plastic bar: PVT wrapped in Al foil then black mylar tape
    G4double bsc_u  = fConfig.backscint_u_cm     * cm;
    G4double bsc_v  = fConfig.backscint_v_cm     * cm;
    G4double bsc_th = fConfig.backscint_thick_cm * cm;
    G4double bsc_gap= fConfig.backscint_gap_cm   * cm;
    // Al envelope half-sizes (Al foil directly on scint surface)
    G4double bscAl_hu = (bsc_u  + 2*tBscAl) / 2;
    G4double bscAl_hv = (bsc_v  + 2*tBscAl) / 2;
    G4double bscAl_hw = (bsc_th + 2*tBscAl) / 2;
    // Tape envelope half-sizes (outer, wraps Al+scint)
    G4double bscTape_hu = (bsc_u  + 2*tBscAl + 2*tTape) / 2;
    G4double bscTape_hv = (bsc_v  + 2*tBscAl + 2*tTape) / 2;
    G4double bscTape_hw = (bsc_th + 2*tBscAl + 2*tTape) / 2;

    // ── World volume ─────────────────────────────────────────
    // Per-axis MM front-face distances (target at mylar-box centre) and the
    // per-arm tangential pinwheel shift.  Arm order: 0=D(+X) 1=B(−X) 2=A(+Z) 3=C(−Z).
    G4double distX     = fConfig.mm_distance_x_cm * cm;
    G4double distZ     = fConfig.mm_distance_z_cm * cm;
    G4double distMax   = std::max(distX, distZ);

    // ── Depth chain (all measured from the MM drift-mylar front face) ────────
    G4double tMM   = tMylar+tAlWin+tKapCath+tCuCath+tDrift+tMesh+tAmp+tResPaste;
    G4double tPCB  = tPCBKap + 4*(tPCBCu+tPCBFR4) + tPCBRoh + tPCBAl;
    G4double mmPcbBack   = tMM + tPCB;                                  // ≈ 3.60 cm
    // SiPM wall: container front measured 11 cm from mylar; scint plane centered.
    G4double sipmFront   = fConfig.sipm_front_from_mylar_cm * cm;       // 11 cm
    G4double sipmContD   = fConfig.sipm_container_depth_cm  * cm;       // 3.3 cm
    G4double sipmScintW  = sipmFront + sipmContD / 2.0;                 // scint depth (center)
    G4double sipmBack    = sipmFront + sipmContD;                       // container back
    // Plastics + LS vessel: per-arm depths (measured 2026-07-17 from the SiPM
    // container back; arm order 0=D 1=B 2=A 3=C).  The LS reference plane is
    // the FLAT slab front face — the front bulge apex sits hCap closer in.
    G4double plasticEnvD = 2 * bscTape_hw;                             // wrapped-bar depth
    G4double plasticWA[4], lsSlabFrontA[4];
    G4double stackDepth = 0.0;                                          // outermost extent
    for (int i = 0; i < 4; ++i) {
        G4double pf = sipmBack + fConfig.gap_sipm_to_plastic_cm[i] * cm;
        plasticWA[i]    = pf + plasticEnvD / 2.0;                       // plastics centre depth
        lsSlabFrontA[i] = sipmBack + fConfig.ls_front_from_sipm_back_cm[i] * cm;
        stackDepth = std::max(stackDepth, lsSlabFrontA[i] + 2*lsTo + hCap);
    }

    G4double lsVExtent   = lsVo + lsFunL + lsNkL + pmtOut;  // vessel+PMT reach along v
    G4double worldHalfXZ = distMax + stackDepth + 5.0*cm;
    G4double worldHalfY  = std::max({sipmBar_hv, lsVExtent, bscTape_hv}) + 5.0*cm;

    auto* worldBox = new G4Box("World", worldHalfXZ, worldHalfY, worldHalfXZ);
    auto* worldLV  = new G4LogicalVolume(worldBox, matAir, "World");
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto* worldPV  = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV,
                                        "World", nullptr, false, 0, true);

    // ── He-3 pressurised capsule at origin ───────────────────
    // Source geometry: STEP file "MASTINU X17 HPRV 00 01 (Cylinder D20 L40 mm)"
    // Polycone axis = local Z; rotateX(-90°) maps local Z → world Y (beam axis).
    //
    // Gas bore:   r=10 mm cylinder (±20 mm) + lower hemispherical end +
    //             conical fill channel up the neck to the Ø1.5 mm valve bore
    // Al vessel:  STEP-derived outer profile, z=−35 mm (tip) to z=+51 mm (valve top)
    //             barrel OD 21.2 mm; shoulder taper into the Ø7 mm neck/valve
    // CFRP wrap:  900 µm over the Al outer surface
    // (gas/Al/CFRP profiles extracted by axial sectioning of the STEP solid;
    //  kept in sync with scripts/plot_geometry.py)

    auto* capRot = new G4RotationMatrix();
    capRot->rotateX(-90.*deg);

    // ── Gas cavity polycone (full STEP-derived interior) ─────
    // He-3 fills the r=10 mm bore (±20 mm), the lower hemispherical end,
    // and the conical fill channel up the neck, necking to the Ø1.5 mm
    // valve bore (r=0.75 mm) at the top.  Profile extracted from the STEP
    // cavity by axial sectioning (see scripts/plot_geometry.py); nests
    // inside the Al vessel with a ≥0.6 mm wall everywhere.
    static const G4int nGas = 23;
    static const G4double zGas[nGas] = {
        -29.500*mm, -28.000*mm, -26.000*mm, -24.000*mm, -22.000*mm,
        -20.000*mm, -15.000*mm,  -5.000*mm,   5.000*mm,  15.000*mm,
         20.000*mm,  22.000*mm,  24.000*mm,  26.000*mm,  28.000*mm,
         30.000*mm,  32.000*mm,  34.000*mm,  36.000*mm,  38.000*mm,
         40.000*mm,  44.000*mm,  50.700*mm
    };
    static const G4double roGas[nGas] = {
          0.001*mm,   6.000*mm,   8.000*mm,   9.165*mm,   9.798*mm,
         10.000*mm,  10.000*mm,  10.000*mm,  10.000*mm,  10.000*mm,
         10.000*mm,   9.798*mm,   9.165*mm,   8.000*mm,   6.299*mm,
          4.842*mm,   3.660*mm,   2.711*mm,   1.967*mm,   1.410*mm,
          1.026*mm,   0.750*mm,   0.750*mm
    };
    static const G4double riGas[nGas] = {
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
    };

    auto* he3Solid = new G4Polycone("He3Gas", 0, 360.*deg, nGas, zGas, riGas, roGas);
    auto* he3LV    = new G4LogicalVolume(he3Solid, GetMat("He3Gas"), "He3Gas");
    he3LV->SetVisAttributes(new G4VisAttributes(G4Color(0.6, 0.9, 1.0, 0.4)));
    he3LV->SetUserLimits(new G4UserLimits(1.0*mm));

    // ── Aluminium vessel polycone (full STEP-derived profile) ─
    // Outer wall sectioned from the STEP at successive z levels; tip at
    // z=−35 mm, barrel OD 21.2 mm (±20 mm), shoulder taper into the
    // Ø7 mm neck/valve (z up to +51 mm).  Matches scripts/plot_geometry.py.
    static const G4int nAl = 29;
    static const G4double zVessel[nAl] = {
        -35.000*mm, -34.000*mm, -33.000*mm, -31.000*mm, -29.000*mm,
        -27.000*mm, -25.000*mm, -23.000*mm, -21.000*mm, -20.000*mm,
        -15.000*mm,  -5.000*mm,   5.000*mm,  15.000*mm,  20.000*mm,
         21.000*mm,  23.000*mm,  25.000*mm,  27.000*mm,  29.000*mm,
         31.000*mm,  33.000*mm,  35.000*mm,  37.000*mm,  39.000*mm,
         40.000*mm,  45.000*mm,  50.000*mm,  51.000*mm
    };
    static const G4double roAl[nAl] = {
          0.000*mm,   3.803*mm,   5.287*mm,   7.206*mm,   8.480*mm,
          9.375*mm,   9.994*mm,  10.386*mm,  10.600*mm,  10.600*mm,
         10.600*mm,  10.600*mm,  10.600*mm,  10.600*mm,  10.600*mm,
         10.600*mm,  10.386*mm,   9.994*mm,   9.375*mm,   8.480*mm,
          7.206*mm,   5.747*mm,   4.708*mm,   4.015*mm,   3.621*mm,
          3.500*mm,   3.500*mm,   3.500*mm,   3.500*mm
    };
    static const G4double riAl[nAl] = {
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
    };

    auto* alSolid = new G4Polycone("He3Cap_Al", 0, 360.*deg, nAl, zVessel, riAl, roAl);
    auto* alLV    = new G4LogicalVolume(alSolid, matAl, "He3Cap_Al");
    alLV->SetVisAttributes(new G4VisAttributes(G4Color(0.7, 0.7, 0.7, 0.8)));

    // ── CFRP wrap polycone (Al outer + 0.9 mm) ───────────────
    static const G4double roCFRP[nAl] = {
          0.000*mm,   4.703*mm,   6.187*mm,   8.106*mm,   9.380*mm,
         10.275*mm,  10.894*mm,  11.286*mm,  11.500*mm,  11.500*mm,
         11.500*mm,  11.500*mm,  11.500*mm,  11.500*mm,  11.500*mm,
         11.500*mm,  11.286*mm,  10.894*mm,  10.275*mm,   9.380*mm,
          8.106*mm,   6.647*mm,   5.608*mm,   4.915*mm,   4.521*mm,
          4.400*mm,   4.400*mm,   4.400*mm,   4.400*mm
    };
    static const G4double riCFRP[nAl] = {
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
    };

    auto* cfrpSolid = new G4Polycone("He3Cap_CFRP", 0, 360.*deg, nAl, zVessel, riCFRP, roCFRP);
    auto* cfrpLV    = new G4LogicalVolume(cfrpSolid, GetMat("CFRP"), "He3Cap_CFRP");
    cfrpLV->SetVisAttributes(new G4VisAttributes(G4Color(0.15, 0.15, 0.15, 0.9)));

    // Place: CFRP in world, Al inside CFRP, gas inside Al
    new G4PVPlacement(capRot, G4ThreeVector(), cfrpLV, "He3Cap_CFRP", worldLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(), alLV,  "He3Cap_Al",   cfrpLV,  false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(), he3LV, "He3Gas",      alLV,    false, 0, true);

    G4cout << "\n=== X17 Full-Experiment Geometry ===" << G4endl;
    G4cout << "  Beam axis    : +Y" << G4endl;
    G4cout << "  He-3 target  : r=10 mm bore, L=40 mm + dome + neck fill channel, 500 bar" << G4endl;
    G4cout << "  Al vessel    : STEP profile, z=-35 to +51 mm (tip to valve)" << G4endl;
    G4cout << "  Gas mixture  : " << fConfig.gas
           << "  (rho=" << matGas->GetDensity()/(mg/cm3) << " mg/cm3)" << G4endl;
    G4cout << "  Stack depth  : " << stackDepth/cm << " cm" << G4endl;
    G4cout << "  LS vessel    : slab " << 2*lsUo/cm << "x" << 2*lsVo/cm << " cm, "
           << 2*lsTo/mm << " mm + bulge domes h=" << hCap/mm << " mm/side; "
           << "funnel " << lsFunL/cm << " cm; neck r=" << lsNkR/mm << " mm x "
           << lsNkL/cm << " cm, PMT in " << pmtIns/cm << " cm"
           << " (fill target " << fConfig.ls_fill_liters << " L)" << G4endl;

    // ── Helper: create a named logical volume (box) ──────────
    auto MakeLV = [&](const std::string& name, G4double hu, G4double hv,
                       G4double hw, G4Material* mat,
                       const G4Color& col, bool solid=false) -> G4LogicalVolume* {
        auto* box = new G4Box(name, hu, hv, hw);
        auto* lv  = new G4LogicalVolume(box, mat, name);
        auto* vis = new G4VisAttributes(col);
        if (solid) vis->SetForceSolid(true);
        lv->SetVisAttributes(vis);
        return lv;
    };

    // ── Create logical volumes for each layer type ───────────
    // (shared across 4 arms via copy number = arm ID)

    // MM layers (38×34 cm face)
    auto* mylarLV   = MakeLV("GasWindow_Mylar",     mmU_hf,mmV_hf,tMylar/2,    matMylar,  G4Color(0.7,0.9,0.7,0.5));
    auto* alWinLV   = MakeLV("GasWindow_Al",         mmU_hf,mmV_hf,tAlWin/2,    matAl,     G4Color(0.7,0.7,0.7,0.8));
    auto* kapCathLV = MakeLV("DriftCathode_Kapton",  mmU_hf,mmV_hf,tKapCath/2,  matKapton, G4Color(0.9,0.7,0.0,0.7));
    auto* cuCathLV  = MakeLV("DriftCathode_Cu",      mmU_hf,mmV_hf,tCuCath/2,   matCu,     G4Color(0.8,0.4,0.1,0.8));
    fDriftGasLV     = MakeLV("DriftGas",             mmU_hf,mmV_hf,tDrift/2,    matGas,    G4Color(0.2,0.5,1.0,0.3));
    auto* meshLV    = MakeLV("Micromesh",             mmU_hf,mmV_hf,tMesh/2,     matSteel,  G4Color(0.5,0.5,0.5,0.9));
    fAmpGasLV       = MakeLV("AmpGas",               mmU_hf,mmV_hf,tAmp/2,      matGas,    G4Color(1.0,0.3,0.3,0.3));
    auto* resLV     = MakeLV("ResistivePaste",        mmU_hf,mmV_hf,tResPaste/2, GetMat("ResistivePaste"),G4Color(0.2,0.2,0.2,0.8));

    // PCB layers (38×34 cm face)
    auto* pcbKapLV  = MakeLV("PCB_Kapton",   mmU_hf,mmV_hf,tPCBKap/2,  matKapton,         G4Color(0.9,0.7,0.0,0.7));
    auto* pcbCu1LV  = MakeLV("PCB_Cu_1",     mmU_hf,mmV_hf,tPCBCu/2,   matCu,             G4Color(0.8,0.4,0.1,0.8));
    auto* pcbFR41LV = MakeLV("PCB_FR4_1",    mmU_hf,mmV_hf,tPCBFR4/2,  GetMat("FR4"),     G4Color(0.2,0.6,0.2,0.8));
    auto* pcbCu2LV  = MakeLV("PCB_Cu_2",     mmU_hf,mmV_hf,tPCBCu/2,   matCu,             G4Color(0.8,0.4,0.1,0.8));
    auto* pcbFR42LV = MakeLV("PCB_FR4_2",    mmU_hf,mmV_hf,tPCBFR4/2,  GetMat("FR4"),     G4Color(0.2,0.6,0.2,0.8));
    auto* pcbCu3LV  = MakeLV("PCB_Cu_3",     mmU_hf,mmV_hf,tPCBCu/2,   matCu,             G4Color(0.8,0.4,0.1,0.8));
    auto* pcbFR43LV = MakeLV("PCB_FR4_3",    mmU_hf,mmV_hf,tPCBFR4/2,  GetMat("FR4"),     G4Color(0.2,0.6,0.2,0.8));
    auto* pcbCu4LV  = MakeLV("PCB_Cu_4",     mmU_hf,mmV_hf,tPCBCu/2,   matCu,             G4Color(0.8,0.4,0.1,0.8));
    auto* pcbFR44LV = MakeLV("PCB_FR4_4",    mmU_hf,mmV_hf,tPCBFR4/2,  GetMat("FR4"),     G4Color(0.2,0.6,0.2,0.8));
    auto* pcbRohLV  = MakeLV("PCB_Rohacell", mmU_hf,mmV_hf,tPCBRoh/2,  GetMat("Rohacell51"),G4Color(0.9,0.9,0.6,0.5));
    auto* pcbAlLV   = MakeLV("PCB_AlFoil",   mmU_hf,mmV_hf,tPCBAl/2,   matAl,             G4Color(0.7,0.7,0.7,0.8));

    // SiPM trigger wall — one plastic scint bar (2.5 cm wide, 50 cm long),
    // placed 16× per arm along u (copyNo = arm).  Scored as "PlasticScint".
    fPlScintLV = MakeLV("PlasticScint", sipmBar_hu, sipmBar_hv, sipmBar_hw,
                         matPlScint, G4Color(0.9,0.9,0.2,0.7));

    // ── LS vessel (STEP "LS X17.step") ───────────────────────
    // CFRP shell = slab + funnel loft (Trd approximation of the STEP B-spline
    // loft) + neck + two ellipsoid bulge domes; the LAB liquid (scored) is the
    // matching interior union placed as a daughter of the shell; the neck bore
    // beyond the liquid holds the half-inserted PMT.  All solids are built in
    // the arm local frame (x=u, y=v with the neck at +v, z=w).  Placements and
    // boolean displacements below use G4Transform3D, which is the DIRECT
    // (active) transform of the solid — unlike the pRot constructor, which
    // rotates the frame.  kToV turns the Trd/Tubs z-axes onto local +v.
    const G4Transform3D kToV = G4RotateX3D(-90.*deg);

    auto* oSlab = new G4Box("LS_oSlab", lsUo, lsVo, lsTo);
    auto* oFun  = new G4Trd("LS_oFunnel", lsUo, lsNkR, lsTo, lsNkR, lsFunL/2);
    auto* oNeck = new G4Tubs("LS_oNeck", 0, lsNkR, lsNkL/2, 0, twopi);
    auto* oCap  = new G4Ellipsoid("LS_oCap", lsUo, lsVo, hCap);
    G4VSolid* lsShellS = new G4UnionSolid("LS_v1", oSlab, oFun,
        G4Translate3D(0, lsVo + lsFunL/2, 0) * kToV);
    lsShellS = new G4UnionSolid("LS_v2", lsShellS, oNeck,
        G4Translate3D(0, lsVo + lsFunL + lsNkL/2, 0) * kToV);
    lsShellS = new G4UnionSolid("LS_v3", lsShellS, oCap, G4Translate3D(0, 0, +lsTo));
    lsShellS = new G4UnionSolid("LS_Vessel", lsShellS, oCap, G4Translate3D(0, 0, -lsTo));

    // Liquid interior — built in the frame of its first solid (the interior
    // slab, whose centre sits at +lsSlabInCen in the vessel frame; the daughter
    // placement below restores that offset).
    G4double stubHalf = (pmtFaceV - stubStartV) / 2;
    auto* iSlab = new G4Box("LS_iSlab", lsUi, lsSlabInHV, lsTi);
    auto* iFun  = new G4Trd("LS_iFunnel", lsUi, iFunU2, lsTi, iFunT2, iFunHL);
    auto* iStub = new G4Tubs("LS_iStub", 0, lsNkRi, stubHalf, 0, twopi);
    auto* iCap  = new G4Ellipsoid("LS_iCap", lsUi, lsSlabInHV, hCap);
    G4VSolid* lsLiqS = new G4UnionSolid("LS_l1", iSlab, iFun,
        G4Translate3D(0, lsVo + iFunHL - lsSlabInCen, 0) * kToV);
    lsLiqS = new G4UnionSolid("LS_l2", lsLiqS, iStub,
        G4Translate3D(0, stubStartV + stubHalf - lsSlabInCen, 0) * kToV);
    lsLiqS = new G4UnionSolid("LS_l3", lsLiqS, iCap, G4Translate3D(0, 0, +lsTi));
    lsLiqS = new G4UnionSolid("LiqScint_1", lsLiqS, iCap, G4Translate3D(0, 0, -lsTi));

    auto* lsShellLV = new G4LogicalVolume(lsShellS, GetMat("CFRP"), "LS_VesselCFRP");
    lsShellLV->SetVisAttributes(new G4VisAttributes(G4Color(0.15,0.15,0.15,0.9)));
    fLS1LV = new G4LogicalVolume(lsLiqS, GetMat("LAB"), "LiqScint_1");
    fLS1LV->SetVisAttributes(new G4VisAttributes(G4Color(0.3,0.8,0.9,0.4)));

    // PMT bore (air) fills the neck from just past the liquid to just short of
    // the neck opening (0.05 mm insets keep daughter surfaces off the shell
    // boundary; the resulting ~0.1 mm CFRP films are negligible).
    G4double boreLen = lsNkL - pmtIns - 0.1*mm;
    auto* boreS  = new G4Tubs("LS_PMTBore", 0, lsNkRi, boreLen/2, 0, twopi);
    auto* boreLV = new G4LogicalVolume(boreS, matAir, "LS_PMTBoreAir");
    boreLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // PMT: borosilicate envelope + vacuum, split at the neck opening into an
    // inserted segment (daughter of the bore) and a protruding segment (world).
    G4Material* matPyrex = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
    G4Material* matVac   = nist->FindOrBuildMaterial("G4_Galactic");
    G4double pmtInsHL = boreLen/2 - 0.05*mm;
    auto* pmtInGlassLV = new G4LogicalVolume(
        new G4Tubs("LS_PMTGlassIn", 0, pmtR, pmtInsHL, 0, twopi),
        matPyrex, "LS_PMTGlassIn");
    auto* pmtInVacLV = new G4LogicalVolume(
        new G4Tubs("LS_PMTVacIn", 0, pmtR - pmtGlass, pmtInsHL - pmtGlass, 0, twopi),
        matVac, "LS_PMTVacIn");
    auto* pmtOutGlassLV = new G4LogicalVolume(
        new G4Tubs("LS_PMTGlassOut", 0, pmtR, pmtOut/2, 0, twopi),
        matPyrex, "LS_PMTGlassOut");
    auto* pmtOutVacLV = new G4LogicalVolume(
        new G4Tubs("LS_PMTVacOut", 0, pmtR - pmtGlass, pmtOut/2 - pmtGlass, 0, twopi),
        matVac, "LS_PMTVacOut");
    for (auto* lv : {pmtInGlassLV, pmtOutGlassLV})
        lv->SetVisAttributes(new G4VisAttributes(G4Color(0.75,0.75,0.78,0.8)));
    for (auto* lv : {pmtInVacLV, pmtOutVacLV})
        lv->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Vessel-internal nesting (shared by all four arms)
    new G4PVPlacement(nullptr, G4ThreeVector(0, lsSlabInCen, 0), fLS1LV,
                      "LiqScint_1_in_vessel", lsShellLV, false, 0, true);
    new G4PVPlacement(G4Translate3D(0, pmtFaceV + 0.05*mm + boreLen/2, 0) * kToV,
                      boreLV, "LS_PMTBore_in_vessel", lsShellLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(), pmtInGlassLV,
                      "LS_PMTGlassIn_in_bore", boreLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(), pmtInVacLV,
                      "LS_PMTVacIn_in_glass", pmtInGlassLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(), pmtOutVacLV,
                      "LS_PMTVacOut_in_glass", pmtOutGlassLV, false, 0, true);

    // Back plastic scint bars: tape (outer, 200µm mylar) → Al foil (20µm) → PVT scint
    // Two LVs per component (L/R) so volume name encodes which bar.
    auto* bscTapeLLV = MakeLV("BackScintTapeL", bscTape_hu,bscTape_hv,bscTape_hw,
                               matBlkMylar, G4Color(0.1,0.1,0.1,0.7));
    auto* bscTapeRLV = MakeLV("BackScintTapeR", bscTape_hu,bscTape_hv,bscTape_hw,
                               matBlkMylar, G4Color(0.1,0.1,0.1,0.7));
    auto* bscAlLLV   = MakeLV("BackScintAlL",   bscAl_hu,  bscAl_hv,  bscAl_hw,
                               matAl,       G4Color(0.7,0.7,0.7,0.6));
    auto* bscAlRLV   = MakeLV("BackScintAlR",   bscAl_hu,  bscAl_hv,  bscAl_hw,
                               matAl,       G4Color(0.7,0.7,0.7,0.6));
    fBackScintLLV    = MakeLV("BackScintL",  bsc_u/2,bsc_v/2,bsc_th/2,
                               matPlScint, G4Color(0.9,0.5,0.1,0.8));
    fBackScintRLV    = MakeLV("BackScintR",  bsc_u/2,bsc_v/2,bsc_th/2,
                               matPlScint, G4Color(0.9,0.5,0.1,0.8));

    // Nesting: tape → Al → PVT (all centred at parent origin)
    new G4PVPlacement(nullptr, G4ThreeVector(), bscAlLLV,    "BackScintAlL_in_tape", bscTapeLLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(), bscAlRLV,    "BackScintAlR_in_tape", bscTapeRLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(), fBackScintLLV,"BackScintL_in_al",    bscAlLLV,   false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(), fBackScintRLV,"BackScintR_in_al",    bscAlRLV,   false, 0, true);

    fDriftGasLV->SetUserLimits(new G4UserLimits(100.*um));
    fAmpGasLV  ->SetUserLimits(new G4UserLimits(100.*um));

    // ── Ordered slab lists (depths accumulate from each group's front) ────────
    struct Slab { G4LogicalVolume* lv; G4double thickness; };
    // MM + PCB stack — placed relative to the arm front face (pinwheel-shifted).
    std::vector<Slab> mmSlabs = {
        // MM
        {mylarLV,   tMylar},
        {alWinLV,   tAlWin},
        {kapCathLV, tKapCath},
        {cuCathLV,  tCuCath},
        {fDriftGasLV, tDrift},
        {meshLV,    tMesh},
        {fAmpGasLV, tAmp},
        {resLV,     tResPaste},
        // PCB
        {pcbKapLV,  tPCBKap},
        {pcbCu1LV,  tPCBCu}, {pcbFR41LV, tPCBFR4},
        {pcbCu2LV,  tPCBCu}, {pcbFR42LV, tPCBFR4},
        {pcbCu3LV,  tPCBCu}, {pcbFR43LV, tPCBFR4},
        {pcbCu4LV,  tPCBCu}, {pcbFR44LV, tPCBFR4},
        {pcbRohLV,  tPCBRoh},
        {pcbAlLV,   tPCBAl},
    };

    // ── SiPM read-out window: which of the 20 bars are instrumented ───────────
    // Bars are centered on the STRUCTURE (u = 0).  The MM sits at local −u
    // (pinwheel shift is along −uHat), so "toward the MM" = −u.  Keep the 16
    // bars whose window centre is offset sipm_readout_shift_bars toward −u:
    // drop 3 bars on the far (+u) side and 1 on the near (−u) side.
    const int    Nb  = fConfig.sipm_n_bars;
    const int    Nr  = fConfig.sipm_n_readout;
    const double bw  = fConfig.sipm_bar_width_cm * cm;
    const double barHalf   = (Nr - 1) / 2.0;
    const double barCenter = (Nb - 1) / 2.0 - fConfig.sipm_readout_shift_bars;
    const int    barLo = static_cast<int>(std::lround(barCenter - barHalf));
    const int    barHi = static_cast<int>(std::lround(barCenter + barHalf));

    // ── Arm placement definitions ─────────────────────────────
    struct ArmDef {
        G4RotationMatrix* rot;
        G4ThreeVector     frontFace;
        G4ThreeVector     uHat;    // local u direction in world
        G4ThreeVector     wHat;    // local w (depth) direction in world
    };

    auto* rot0 = new G4RotationMatrix(); rot0->rotateY( 90.*deg);
    auto* rot1 = new G4RotationMatrix(); rot1->rotateY(-90.*deg);
    auto* rot3 = new G4RotationMatrix(); rot3->rotateY(180.*deg);

    ArmDef armDefs[4] = {
        {rot0, {distX,0,0},  {0,0,-1}, {1,0,0}},   // Arm0 +X = D
        {rot1, {-distX,0,0}, {0,0, 1}, {-1,0,0}},  // Arm1 −X = B
        {nullptr, {0,0,distZ}, {1,0,0}, {0,0,1}},  // Arm2 +Z = A
        {rot3, {0,0,-distZ}, {-1,0,0}, {0,0,-1}},  // Arm3 −Z = C
    };

    // ── Place arms ────────────────────────────────────────────
    for (int arm = 0; arm < 4; ++arm) {
        const auto& ad = armDefs[arm];

        // Pinwheel: slide the whole arm tangentially (⟂ its normal) along
        // −uHat by the measured per-MM amount.  This translates the arm's
        // local coordinate origin too, so hit u-coordinates stay centred.
        G4ThreeVector pinShift = -ad.uHat * (fConfig.mm_pinwheel_shift_cm[arm] * cm);
        G4ThreeVector armFront = ad.frontFace + pinShift;

        // The SiPM wall is centred on the mechanical STRUCTURE (un-shifted),
        // while MM / plastics / LS are centred on the pinwheel-shifted MM.
        G4ThreeVector structFront = ad.frontFace;

        // Store arm axes for SteppingAction coordinate transforms (MM frame)
        fArmAxes[arm].frontFace = armFront;
        fArmAxes[arm].uHat      = ad.uHat;
        fArmAxes[arm].vHat      = G4ThreeVector(0, 1, 0);
        fArmAxes[arm].wHat      = ad.wHat;

        // Helper: place lv at local (u, 0, w) relative to `base` (arm front).
        auto place = [&](G4LogicalVolume* lv, const G4ThreeVector& base,
                          G4double u, G4double w, const std::string& tag) {
            G4ThreeVector localPos(u, 0, w);
            G4ThreeVector worldPos = base +
                (ad.rot ? (*ad.rot)*localPos : localPos);
            new G4PVPlacement(ad.rot, worldPos, lv,
                               "Arm" + std::to_string(arm) + "_" + tag,
                               worldLV, false, arm, true);
        };

        // 1) MM + PCB stack — relative to the (shifted) arm front face.
        G4double zLocal = 0.0;
        for (const auto& s : mmSlabs) {
            zLocal += s.thickness;
            place(s.lv, armFront, 0.0, zLocal - s.thickness / 2.0, s.lv->GetName());
        }

        // 2) SiPM wall — 16 instrumented bars, centred on the STRUCTURE.
        for (int i = barLo; i <= barHi; ++i) {
            G4double u_i = bw * (i - (Nb - 1) / 2.0);   // bar centre from structure centre
            place(fPlScintLV, structFront, u_i, sipmScintW,
                  "PlasticScint_bar" + std::to_string(i));
        }

        // 3) Plastics — two wrapped bars side-by-side, centred on the MM;
        //    per-arm front distance (measured 2026-07-17).
        G4double uOff = bscTape_hu + bsc_gap / 2.0;
        place(bscTapeLLV, armFront, -uOff, plasticWA[arm], "BackTapeL");
        place(bscTapeRLV, armFront, +uOff, plasticWA[arm], "BackTapeR");

        // 4) LS vessel — surveyed 2026-07-17: flat slab front face at the
        //    measured per-arm depth, slab centre at the surveyed height
        //    (ls_offset_v), slab centred on the MM in u (ASSUMED, not
        //    measured).  ls_rot_deg rotates the vessel about its depth axis:
        //    0 = vertical, neck/PMT up (B, C); −90 = horizontal, neck/PMT
        //    along +u (A, D).  Placed with G4Transform3D (the direct/active-
        //    rotation constructor), matching the active use of ad.rot in the
        //    position math above.
        G4RotationMatrix lsRz; lsRz.rotateZ(fConfig.ls_rot_deg[arm] * deg);
        G4RotationMatrix lsArmR = (ad.rot ? *ad.rot : G4RotationMatrix()) * lsRz;
        G4double lsSlabCenW = lsSlabFrontA[arm] + lsTo;
        G4ThreeVector lsLocal(0, fConfig.ls_offset_v_cm[arm] * cm, lsSlabCenW);
        G4ThreeVector lsOrigin = armFront + (ad.rot ? (*ad.rot)*lsLocal : lsLocal);
        new G4PVPlacement(G4Transform3D(lsArmR, lsOrigin), lsShellLV,
                          "Arm" + std::to_string(arm) + "_LS_Vessel",
                          worldLV, false, arm, true);
        if (pmtOut > 0.1*mm) {   // protruding PMT segment beyond the neck opening
            G4RotationMatrix rxm90; rxm90.rotateX(-90.*deg);   // solid z → vessel +y
            G4ThreeVector pmtPos = lsOrigin +
                lsArmR * G4ThreeVector(0, lsVo + lsFunL + lsNkL + pmtOut/2, 0);
            new G4PVPlacement(G4Transform3D(lsArmR * rxm90, pmtPos), pmtOutGlassLV,
                              "Arm" + std::to_string(arm) + "_LS_PMT",
                              worldLV, false, arm, true);
        }

        G4cout << "  Arm " << arm << " front face at " << armFront/cm
               << " cm  (pinwheel " << fConfig.mm_pinwheel_shift_cm[arm]
               << " cm); SiPM bars " << barLo << "-" << barHi << G4endl;
    }

    G4cout << "=====================================" << G4endl;
    return worldPV;
}

// ─────────────────────────────────────────────────────────────────────────────
void DetectorConstruction::ConstructSDandField() {
    auto RegisterSD = [&](G4LogicalVolume* lv, const std::string& sdName) {
        if (!lv) return;
        auto* sd = new SensitiveDetector(sdName, sdName + "Hits",
                                          lv->GetName(), fConfig);
        G4SDManager::GetSDMpointer()->AddNewDetector(sd);
        SetSensitiveDetector(lv, sd);
    };
    RegisterSD(fDriftGasLV,   "DriftGasSD");
    RegisterSD(fAmpGasLV,     "AmpGasSD");
    RegisterSD(fPlScintLV,    "PlasticScintSD");
    RegisterSD(fLS1LV,        "LiqScint1SD");
    RegisterSD(fBackScintLLV, "BackScintLSD");
    RegisterSD(fBackScintRLV, "BackScintRSD");
}
