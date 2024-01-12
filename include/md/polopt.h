// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

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
MDQC_EXTERN std::vector<real> copt;
MDQC_EXTERN std::vector<real> copm;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> uopt;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> uoptp;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> uopts;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> uoptps;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> fopt;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> foptp;
}
