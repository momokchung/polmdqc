/////////////////////////////////////////////////////////
//                                                     //
//  bound.h  --  periodic boundary condition controls  //
//                                                     //
/////////////////////////////////////////////////////////

// "sizes" sets values for array dimensions used throughout
// the software; these parameters fix the size of the largest
// systems that can be handled

// polycut       cutoff distance for infinite polymer nonbonds
// polycut2      square of infinite polymer nonbond cutoff
// use_bounds    flag to use periodic boundary conditions
// use_replica   flag to use replicates for periodic system
// use_polymer   flag to mark presence of infinite polymer


#pragma once
#include "macro.h"
#include <string>

QCMD_EXTERN double polycut;
QCMD_EXTERN double polycut2;
QCMD_EXTERN bool use_bounds;
QCMD_EXTERN bool use_replica;
QCMD_EXTERN bool use_polymer;
