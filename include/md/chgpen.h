// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  chgpen  --  charge penetration in current structure  //
//                                                       //
///////////////////////////////////////////////////////////

// ncp       total number of charge penetration sites in system
// pcore     number of core electrons assigned to each atom
// pval      number of valence electrons assigned to each atom
// pval0     original number of valence electrons for charge flux
// palpha    charge penetration damping value at each atom

MDQC_EXTERN int ncp;
MDQC_EXTERN std::vector<double> pcore;
MDQC_EXTERN std::vector<double> pval;
MDQC_EXTERN std::vector<double> pval0;
MDQC_EXTERN std::vector<double> palpha;
}
