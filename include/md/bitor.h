// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  bitor  --  bitorsions in the current structure  //
//                                                  //
//////////////////////////////////////////////////////

// nbitor  total number of bitorsions in the system
// ibitor  numbers of the atoms in each bitorsion

MDQC_EXTERN int nbitor;
MDQC_EXTERN MDQCArray2D<int,5> ibitor;
}
