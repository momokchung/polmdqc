// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  rigid  --  rigid body coordinates for atom groups  //
//                                                     //
/////////////////////////////////////////////////////////

// xrb         rigid body reference x-coordinate for each atom
// yrb         rigid body reference y-coordinate for each atom
// zrb         rigid body reference z-coordinate for each atom
// rbc         current rigid body coordinates for each group
// use_rigid   flag to mark use of rigid body coordinate system 

MDQC_EXTERN std::vector<double> xrb;
MDQC_EXTERN std::vector<double> yrb;
MDQC_EXTERN std::vector<double> zrb;
MDQC_EXTERN std::vector<std::vector<double>> rbc;
MDQC_EXTERN bool use_rigid;
}
