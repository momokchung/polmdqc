// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  katoms  --  atom definition forcefield parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// atmcls     atom class number for each of the atom types
// atmnum     atomic number for each of the atom types
// ligand     number of atoms to be attached to each atom type
// weight     average atomic mass of each atom type
// symbol     modified atomic symbol for each atom type
// describe   string identifying each of the atom types

MDQC_EXTERN MDQCArray<int> atmcls;
MDQC_EXTERN MDQCArray<int> atmnum;
MDQC_EXTERN MDQCArray<int> ligand;
MDQC_EXTERN MDQCArray<real> weight;
MDQC_EXTERN MDQCArray<std::string> symbol;
MDQC_EXTERN MDQCArray<std::string> describe;
}
