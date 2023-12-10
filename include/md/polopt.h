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
MDQC_EXTERN std::vector<double> copt;
MDQC_EXTERN std::vector<double> copm;
MDQC_EXTERN std::vector<std::vector<std::vector<double>>> uopt;
MDQC_EXTERN std::vector<std::vector<std::vector<double>>> uoptp;
MDQC_EXTERN std::vector<std::vector<std::vector<double>>> uopts;
MDQC_EXTERN std::vector<std::vector<std::vector<double>>> uoptps;
MDQC_EXTERN std::vector<std::vector<std::vector<double>>> fopt;
MDQC_EXTERN std::vector<std::vector<std::vector<double>>> foptp;
}
