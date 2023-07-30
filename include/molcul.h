////////////////////////////////////////////////////////////
//                                                        //
//  molcul.h  --  individual molecules in current system  //
//                                                        //
////////////////////////////////////////////////////////////

// nmol      total number of separate molecules in the system
// imol      first and last atom of each molecule in the list
// kmol      contiguous list of the atoms in each molecule
// molcule   number of the molecule to which each atom belongs
// totmass   total weight of all the molecules in the system
// molmass   molecular weight for each molecule in the system


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nmol;
QCMD_EXTERN std::vector<std::vector<int>> imol;
QCMD_EXTERN std::vector<int> kmol;
QCMD_EXTERN std::vector<int> molcule;
QCMD_EXTERN double totmass;
QCMD_EXTERN std::vector<double> molmass;
