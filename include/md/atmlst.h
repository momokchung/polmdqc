// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  atmlst  --  bond and angle local geometry indices  //
//                                                     //
/////////////////////////////////////////////////////////

// bndlist   numbers of the bonds involving each atom
// anglist   numbers of the angles centered on each atom
// balist    numbers of the bonds comprising each angle

MDQC_EXTERN std::vector<std::vector<int>> bndlist;
MDQC_EXTERN std::vector<std::vector<int>> anglist;
MDQC_EXTERN std::vector<std::vector<int>> balist;
}
