/////////////////////////////////////////////////////////
//                                                     //
//  tarray.h  --  store dipole-dipole matrix elements  //
//                                                     //
/////////////////////////////////////////////////////////

// ntpair     number of stored dipole-dipole matrix elements
// tindex     index into stored dipole-dipole matrix values
// tdipdip    stored dipole-dipole matrix element values


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int ntpair;
QCMD_EXTERN std::vector<std::vector<int>> tindex;
QCMD_EXTERN std::vector<std::vector<double>> tdipdip;
