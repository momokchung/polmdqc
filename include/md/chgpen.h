// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

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
MDQC_EXTERN MDQCArray<real> pcore;
MDQC_EXTERN MDQCArray<real> pval;
MDQC_EXTERN MDQCArray<real> pval0;
MDQC_EXTERN MDQCArray<real> palpha;
}
