///////////////////////////////////////////////////////////
//                                                       //
//  usage.h  --  atoms active during energy computation  //
//                                                       //
///////////////////////////////////////////////////////////

// nuse   total number of active atoms in energy calculation
// iuse   numbers of the atoms active in energy calculation
// use    true if an atom is active, false if inactive


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nuse;
QCMD_EXTERN std::vector<int> iuse;
QCMD_EXTERN std::vector<bool> use;
