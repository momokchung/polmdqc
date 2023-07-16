///////////////////////////////////////////////////////////
//                                                       //
//  katoms.h  --  atom definition forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// atmcls     atom class number for each of the atom types
// atmnum     atomic number for each of the atom types
// ligand     number of atoms to be attached to each atom type
// weight     average atomic mass of each atom type
// symbol     modified atomic symbol for each atom type
// describe   string identifying each of the atom types


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN std::vector<int> atmcls;
QCMD_EXTERN std::vector<int> atmnum;
QCMD_EXTERN std::vector<int> ligand;
QCMD_EXTERN std::vector<double> weight;
QCMD_EXTERN std::vector<std::string> symbol;
QCMD_EXTERN std::vector<std::string> describe;
