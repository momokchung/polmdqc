// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  kopdst  --  out-of-plane distance forcefield params  //
//                                                       //
///////////////////////////////////////////////////////////

// maxnopd   maximum number of out-of-plane distance entries
// opds      force constant parameters for out-of-plane distance
// kopd      string of atom classes for out-of-plane distance

MDQC_EXTERN int maxnopd;
MDQC_EXTERN MDQCArray<real> opds;
MDQC_EXTERN MDQCArray<std::string> kopd;
}
