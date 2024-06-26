// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  khbond  --  H-bonding term forcefield parameters  //
//                                                    //
////////////////////////////////////////////////////////

// maxnhb   maximum number of hydrogen bonding pair entries
// radhb    radius parameter for hydrogen bonding pairs
// epshb    well depth parameter for hydrogen bonding pairs
// khb      string of atom types for hydrogen bonding pairs

MDQC_EXTERN int maxnhb;
MDQC_EXTERN MDQCArray<real> radhb;
MDQC_EXTERN MDQCArray<real> epshb;
MDQC_EXTERN MDQCArray<std::string> khb;
}
