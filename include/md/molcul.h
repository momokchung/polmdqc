// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  molcul  --  individual molecules in current system  //
//                                                      //
//////////////////////////////////////////////////////////

// nmol      total number of separate molecules in the system
// imol      first and last atom of each molecule in the list
// kmol      contiguous list of the atoms in each molecule
// molcule   number of the molecule to which each atom belongs
// totmass   total weight of all the molecules in the system
// molmass   molecular weight for each molecule in the system

MDQC_EXTERN int nmol;
MDQC_EXTERN std::vector<std::vector<int>> imol;
MDQC_EXTERN std::vector<int> kmol;
MDQC_EXTERN std::vector<int> molcule;
MDQC_EXTERN real totmass;
MDQC_EXTERN std::vector<real> molmass;
}
