// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  kpolpr  --  special Thole forcefield parameters  //
//                                                   //
///////////////////////////////////////////////////////

// maxnpp   maximum number of special pair polarization entries
// thlpr    Thole damping values for special polarization pairs
// thdpr    Thole direct damping for special polarization pairs
// kppr     string of atom types for special polarization pairs

MDQC_EXTERN int maxnpp;
MDQC_EXTERN MDQCArray<real> thlpr;
MDQC_EXTERN MDQCArray<real> thdpr;
MDQC_EXTERN MDQCArray<std::string> kppr;
}
