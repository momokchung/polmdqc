///////////////////////////////////////////////////////
//                                                   //
//  atoms.h  --  number, position and type of atoms  //
//                                                   //
///////////////////////////////////////////////////////

// n       total number of atoms in the current system
// type    atom type number for each atom in the system
// x       current x-coordinate for each atom in the system
// y       current y-coordinate for each atom in the system
// z       current z-coordinate for each atom in the system


#pragma once
#include "macro.h"
#include "sizes.h"
#include <string>

QCMD_EXTERN int n;
QCMD_EXTERN int type[maxatm];
QCMD_EXTERN double x[maxatm];
QCMD_EXTERN double y[maxatm];
QCMD_EXTERN double z[maxatm];
