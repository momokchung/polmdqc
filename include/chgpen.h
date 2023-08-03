/////////////////////////////////////////////////////////////
//                                                         //
//  chgpen.h  --  charge penetration in current structure  //
//                                                         //
/////////////////////////////////////////////////////////////

// ncp       total number of charge penetration sites in system
// pcore     number of core electrons assigned to each atom
// pval      number of valence electrons assigned to each atom
// pval0     original number of valence electrons for charge flux
// palpha    charge penetration damping value at each atom


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int ncp;
QCMD_EXTERN std::vector<double> pcore;
QCMD_EXTERN std::vector<double> pval;
QCMD_EXTERN std::vector<double> pval0;
QCMD_EXTERN std::vector<double> palpha;
