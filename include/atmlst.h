///////////////////////////////////////////////////////////
//                                                       //
//  atmlst.h  --  bond and angle local geometry indices  //
//                                                       //
///////////////////////////////////////////////////////////

// bndlist   numbers of the bonds involving each atom
// anglist   numbers of the angles centered on each atom
// balist    numbers of the bonds comprising each angle


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<std::vector<int>> bndlist;
QCMD_EXTERN std::vector<std::vector<int>> anglist;
QCMD_EXTERN std::vector<std::vector<int>> balist;
