// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  bndpot  --  bond stretch functional form details  //
//                                                    //
////////////////////////////////////////////////////////

// cbnd      cubic coefficient in bond stretch potential
// qbnd      quartic coefficient in bond stretch potential
// bndunit   convert bond stretch energy to kcal/mole
// bndtyp    type of bond stretch potential energy function

MDQC_EXTERN double cbnd;
MDQC_EXTERN double qbnd;
MDQC_EXTERN double bndunit;
MDQC_EXTERN std::string bndtyp;
}
