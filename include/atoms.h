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
#include "sizes.h"
#include <string>

extern int n;
extern int type[maxatm];
extern double x[maxatm];
extern double y[maxatm];
extern double z[maxatm];
