// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

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

constexpr real epso = 0.1100;
constexpr real epsh = 0.0135;
constexpr real rmino = 1.7025;
constexpr real rminh = 1.3275;
constexpr real awater = 0.033428;
constexpr real slevy = 1.0;
constexpr real shctd = 0.75;
constexpr real cavoff = 0.0;
constexpr real dspoff = 1.056;
MDQC_EXTERN real solvprs,surften;
MDQC_EXTERN real spcut,spoff;
MDQC_EXTERN real stcut,stoff;
MDQC_EXTERN MDQCArray<real> radcav;
MDQC_EXTERN MDQCArray<real> raddsp;
MDQC_EXTERN MDQCArray<real> epsdsp;
MDQC_EXTERN MDQCArray<real> cdsp;
}
