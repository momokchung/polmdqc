//////////////////////////////////////////////////////
//                                                  //
//  couple.h  --  atom neighbor connectivity lists  //
//                                                  //
//////////////////////////////////////////////////////

// n12      number of atoms directly bonded to each atom
// n13      number of atoms in a 1-3 relation to each atom
// n14      number of atoms in a 1-4 relation to each atom
// n15      number of atoms in a 1-5 relation to each atom
// i12      atom numbers of atoms 1-2 connected to each atom
// i13      atom numbers of atoms 1-3 connected to each atom
// i14      atom numbers of atoms 1-4 connected to each atom
// i15      atom numbers of atoms 1-5 connected to each atom


#pragma once
#include "macro.h"
#include "sizes.h"
#include <vector>

QCMD_EXTERN int n12[maxatm];
QCMD_EXTERN std::vector<int> n13;
QCMD_EXTERN std::vector<int> n14;
QCMD_EXTERN std::vector<int> n15;
QCMD_EXTERN int i12[maxatm][maxval];
QCMD_EXTERN std::vector<std::vector<int>> i13;
QCMD_EXTERN std::vector<std::vector<int>> i14;
QCMD_EXTERN std::vector<std::vector<int>> i15;
