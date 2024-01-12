// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  kpitor  --  pi-system torsion forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// maxnpt   maximum number of pi-system torsion parameter entries
// ptcon    force constant parameters for pi-system torsions
// kpt      string of atom classes for pi-system torsion terms

MDQC_EXTERN int maxnpt;
MDQC_EXTERN std::vector<real> ptcon;
MDQC_EXTERN std::vector<std::string> kpt;
}
