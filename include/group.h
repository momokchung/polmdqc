////////////////////////////////////////////////////////////
//                                                        //
//  group.h  --  partitioning of system into atom groups  //
//                                                        //
////////////////////////////////////////////////////////////

// ngrp        total number of atom groups in the system
// kgrp        contiguous list of the atoms in each group
// grplist     number of the group to which each atom belongs
// igrp        first and last atom of each group in the list
// grpmass     total mass of all the atoms in each group
// wgrp        weight for each set of group-group interactions (startIndex = 1)
// use_group   flag to use partitioning of system into groups
// use_intra   flag to include only intragroup interactions
// use_inter   flag to include only intergroup interactions


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int ngrp;
QCMD_EXTERN std::vector<int> kgrp;
QCMD_EXTERN std::vector<int> grplist;
QCMD_EXTERN std::vector<std::vector<int>> igrp;
QCMD_EXTERN std::vector<double> grpmass;
QCMD_EXTERN std::vector<std::vector<double>> wgrp;
QCMD_EXTERN bool use_group;
QCMD_EXTERN bool use_intra;
QCMD_EXTERN bool use_inter;
