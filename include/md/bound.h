// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  bound  --  periodic boundary condition controls  //
//                                                   //
///////////////////////////////////////////////////////

// "sizes" sets values for array dimensions used throughout
// the software; these parameters fix the size of the largest
// systems that can be handled

// polycut       cutoff distance for infinite polymer nonbonds
// polycut2      square of infinite polymer nonbond cutoff
// use_bounds    flag to use periodic boundary conditions
// use_replica   flag to use replicates for periodic system
// use_polymer   flag to mark presence of infinite polymer

MDQC_EXTERN double polycut;
MDQC_EXTERN double polycut2;
MDQC_EXTERN bool use_bounds;
MDQC_EXTERN bool use_replica;
MDQC_EXTERN bool use_polymer;
}
