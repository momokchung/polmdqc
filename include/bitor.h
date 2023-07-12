////////////////////////////////////////////////////////
//                                                    //
//  bitor.h  --  bitorsions in the current structure  //
//                                                    //
////////////////////////////////////////////////////////

// nbitor  total number of bitorsions in the system
// ibitor  numbers of the atoms in each bitorsion


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nbitor;
QCMD_EXTERN std::vector<std::vector<int>> ibitor;
