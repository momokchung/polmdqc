// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  kchrge  --  partial charge forcefield parameters  //
//                                                    //
////////////////////////////////////////////////////////

// chg   partial charge parameters for each atom type

MDQC_EXTERN std::vector<real> chg;
}
