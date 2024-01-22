// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

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

MDQC_EXTERN MDQCArray<real> xrb;
MDQC_EXTERN MDQCArray<real> yrb;
MDQC_EXTERN MDQCArray<real> zrb;
MDQC_EXTERN MDQCArray2D<real,6> rbc;
MDQC_EXTERN bool use_rigid;
}
