// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "sizes.h"

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  group  --  partitioning of system into atom groups  //
//                                                      //
//////////////////////////////////////////////////////////

// ngrp        total number of atom groups in the system
// kgrp        contiguous list of the atoms in each group
// grplist     number of the group to which each atom belongs
// igrp        first and last atom of each group in the list
// grpmass     total mass of all the atoms in each group
// wgrp        weight for each set of group-group interactions (startIndex = 1)
// use_group   flag to use partitioning of system into groups
// use_intra   flag to include only intragroup interactions
// use_inter   flag to include only intergroup interactions

MDQC_EXTERN int ngrp;
MDQC_EXTERN MDQCArray<int> kgrp;
MDQC_EXTERN MDQCArray<int> grplist;
MDQC_EXTERN MDQCArray2D<int,2> igrp;
MDQC_EXTERN MDQCArray<real> grpmass;
MDQC_EXTERN MDQCArray2D<real,maxgrp+1> wgrp;
MDQC_EXTERN bool use_group;
MDQC_EXTERN bool use_intra;
MDQC_EXTERN bool use_inter;
}
