// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  kurybr  --  Urey-Bradley term forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// maxnu   maximum number of Urey-Bradley parameter entries
// ucon    force constant parameters for Urey-Bradley terms
// dst13   ideal 1-3 distance parameters for Urey-Bradley terms
// ku      string of atom classes for Urey-Bradley terms

MDQC_EXTERN int maxnu;
MDQC_EXTERN MDQCArray<real> ucon;
MDQC_EXTERN MDQCArray<real> dst13;
MDQC_EXTERN MDQCArray<std::string> ku;
}
