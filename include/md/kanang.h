// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  kanang  --  angle-angle term forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// anan   angle-angle cross term parameters for each atom class

MDQC_EXTERN std::vector<std::vector<real>> anan;
}
