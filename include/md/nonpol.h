// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  nonpol  --  nonpolar cavity & dispersion parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// epso      water oxygen eps for implicit dispersion term
// epsh      water hydrogen eps for implicit dispersion term
// rmino     water oxygen Rmin for implicit dispersion term
// rminh     water hydrogen Rmin for implicit dispersion term
// awater    water number density at standard temp & pressure
// slevy     enthalpy-to-free energy scale factor for dispersion
// shctd     HCT overlap scale factor for the dispersion integral
// cavoff    radius offset for use in computing cavitation energy
// dspoff    radius offset for the start of dispersion integral
// solvprs   limiting microscopic solvent pressure value
// surften   limiting macroscopic surface tension value
// spcut     starting radius for solvent pressure tapering
// spoff     cutoff radius for solvent pressure tapering
// stcut     starting radius for surface tension tapering
// stoff     cutoff radius for surface tension tapering
// radcav    atomic radius of each atom for cavitation energy
// raddsp    atomic radius of each atom for dispersion energy
// epsdsp    vdw well depth of each atom for dispersion energy
// cdsp      maximum dispersion energy for each atom

constexpr double epso = 0.1100;
constexpr double epsh = 0.0135;
constexpr double rmino = 1.7025;
constexpr double rminh = 1.3275;
constexpr double awater = 0.033428;
constexpr double slevy = 1.0;
constexpr double shctd = 0.75;
constexpr double cavoff = 0.0;
constexpr double dspoff = 1.056;
MDQC_EXTERN double solvprs,surften;
MDQC_EXTERN double spcut,spoff;
MDQC_EXTERN double stcut,stoff;
MDQC_EXTERN std::vector<double> radcav;
MDQC_EXTERN std::vector<double> raddsp;
MDQC_EXTERN std::vector<double> epsdsp;
MDQC_EXTERN std::vector<double> cdsp;
}
