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

MDQC_EXTERN real aewald;
MDQC_EXTERN real aeewald;
MDQC_EXTERN real apewald;
MDQC_EXTERN real adewald;
MDQC_EXTERN std::string boundary;
}
