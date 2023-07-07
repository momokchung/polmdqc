/////////////////////////////////////////////////////////
//                                                     //
//  zclose.h  --  Z-matrix ring openings and closures  //
//                                                     //
/////////////////////////////////////////////////////////

// nadd   number of added bonds between Z-matrix atoms
// ndel   number of bonds between Z-matrix bonds to delete
// iadd   numbers of the atom pairs defining added bonds
// idel   numbers of the atom pairs defining deleted bonds


#pragma once
#include "macro.h"
#include "sizes.h"

QCMD_EXTERN int nadd,ndel;
QCMD_EXTERN int iadd[maxatm][2];
QCMD_EXTERN int idel[maxatm][2];
