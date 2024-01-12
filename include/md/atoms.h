// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "sizes.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  atoms  --  number, position and type of atoms  //
//                                                 //
/////////////////////////////////////////////////////

// n       total number of atoms in the current system
// type    atom type number for each atom in the system
// x       current x-coordinate for each atom in the system
// y       current y-coordinate for each atom in the system
// z       current z-coordinate for each atom in the system

MDQC_EXTERN int n;
MDQC_EXTERN int type[maxatm];
MDQC_EXTERN real x[maxatm];
MDQC_EXTERN real y[maxatm];
MDQC_EXTERN real z[maxatm];
}
