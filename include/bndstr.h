/////////////////////////////////////////////////////////////
//                                                         //
//  bndstr.h  --  bond stretches in the current structure  //
//                                                         //
/////////////////////////////////////////////////////////////

// nbond   total number of bond stretches in the system
// ibnd    numbers of the atoms in each bond stretch
// bk      bond stretch force constants (kcal/mole/Ang**2)
// bl      ideal bond length values in Angstroms


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nbond;
QCMD_EXTERN std::vector<std::vector<int>> ibnd;
QCMD_EXTERN std::vector<double> bk;
QCMD_EXTERN std::vector<double> bl;
