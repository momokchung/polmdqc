///////////////////////////////////////////////////////////
//                                                       //
//  rigid.h  --  rigid body coordinates for atom groups  //
//                                                       //
///////////////////////////////////////////////////////////

// xrb         rigid body reference x-coordinate for each atom
// yrb         rigid body reference y-coordinate for each atom
// zrb         rigid body reference z-coordinate for each atom
// rbc         current rigid body coordinates for each group
// use_rigid   flag to mark use of rigid body coordinate system 


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<double> xrb;
QCMD_EXTERN std::vector<double> yrb;
QCMD_EXTERN std::vector<double> zrb;
QCMD_EXTERN std::vector<std::vector<double>> rbc;
QCMD_EXTERN bool use_rigid;
