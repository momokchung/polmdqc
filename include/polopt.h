///////////////////////////////////////////////////////////
//                                                       //
//  polopt.h  --  induced dipoles for OPT extrapolation  //
//                                                       //
///////////////////////////////////////////////////////////

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


#pragma once
#include "macro.h"
#include <vector>

constexpr int maxopt = 6;
QCMD_EXTERN int optorder;
QCMD_EXTERN int optlevel;
QCMD_EXTERN std::vector<double> copt;
QCMD_EXTERN std::vector<double> copm;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> uopt;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> uoptp;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> uopts;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> uoptps;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> fopt;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> foptp;
