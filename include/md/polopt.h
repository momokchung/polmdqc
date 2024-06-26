// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  polopt  --  induced dipoles for OPT extrapolation  //
//                                                     //
/////////////////////////////////////////////////////////

// maxopt    maximum order for OPT induced dipole extrapolation
// optorder  highest coefficient order for OPT dipole extrapolation
// optlevel  current OPT order for reciprocal potential and field
// copt      coefficients for OPT total induced dipole moments
// copm      coefficients for OPT incremental induced dipole moments
// uopt      OPT induced dipole components at each multipole site
// uoptp     OPT induced dipoles in field used for energy terms
// uopts     OPT GK or PB induced dipoles at each multipole site
// uoptps    OPT induced dipoles in field used for GK or PB energy
// fopt      OPT fractional reciprocal potentials at multipole sites
// foptp     OPT fractional reciprocal potentials for energy terms

constexpr int maxopt = 6;
MDQC_EXTERN int optorder;
MDQC_EXTERN int optlevel;
MDQC_EXTERN MDQCArray<real> copt;
MDQC_EXTERN MDQCArray<real> copm;
MDQC_EXTERN MDQCArray2D<real,3> uopt;
MDQC_EXTERN MDQCArray2D<real,3> uoptp;
MDQC_EXTERN MDQCArray2D<real,3> uopts;
MDQC_EXTERN MDQCArray2D<real,3> uoptps;
MDQC_EXTERN MDQCArray2D<real,10> fopt;
MDQC_EXTERN MDQCArray2D<real,10> foptp;
}
