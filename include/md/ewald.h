// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  ewald  --  Ewald summation parameters and options  //
//                                                     //
/////////////////////////////////////////////////////////

// aewald     current value of Ewald convergence coefficient
// aeewald    Ewald convergence coefficient for electrostatics
// apewald    Ewald convergence coefficient for polarization
// adewald    Ewald convergence coefficient for dispersion
// boundary   Ewald boundary condition; none, tinfoil or vacuum

MDQC_EXTERN double aewald;
MDQC_EXTERN double aeewald;
MDQC_EXTERN double apewald;
MDQC_EXTERN double adewald;
MDQC_EXTERN std::string boundary;
}
