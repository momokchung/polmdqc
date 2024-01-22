// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  usage  --  atoms active during energy computation  //
//                                                     //
/////////////////////////////////////////////////////////

// nuse   total number of active atoms in energy calculation
// iuse   numbers of the atoms active in energy calculation
// use    true if an atom is active, false if inactive (starting index=1)

MDQC_EXTERN int nuse;
MDQC_EXTERN MDQCArray<int> iuse;
MDQC_EXTERN MDQCArray<bool> use;
}
