// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "sizes.h"

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  zclose  --  Z-matrix ring openings and closures  //
//                                                   //
///////////////////////////////////////////////////////

// nadd   number of added bonds between Z-matrix atoms
// ndel   number of bonds between Z-matrix bonds to delete
// iadd   numbers of the atom pairs defining added bonds
// idel   numbers of the atom pairs defining deleted bonds

MDQC_EXTERN int nadd,ndel;
MDQC_EXTERN int iadd[maxatm][2];
MDQC_EXTERN int idel[maxatm][2];
}
