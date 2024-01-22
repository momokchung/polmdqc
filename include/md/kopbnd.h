// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  kopbnd  --  out-of-plane bend forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// maxnopb   maximum number of out-of-plane bending entries
// opbn      force constant parameters for out-of-plane bending
// kopb      string of atom classes for out-of-plane bending

MDQC_EXTERN int maxnopb;
MDQC_EXTERN MDQCArray<real> opbn;
MDQC_EXTERN MDQCArray<std::string> kopb;
}
