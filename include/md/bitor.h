// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  bitor  --  bitorsions in the current structure  //
//                                                  //
//////////////////////////////////////////////////////

// nbitor  total number of bitorsions in the system
// ibitor  numbers of the atoms in each bitorsion

MDQC_EXTERN int nbitor;
MDQC_EXTERN std::vector<std::vector<int>> ibitor;
}
