// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  tarray  --  store dipole-dipole matrix elements  //
//                                                   //
///////////////////////////////////////////////////////

// ntpair     number of stored dipole-dipole matrix elements
// tindex     index into stored dipole-dipole matrix values
// tdipdip    stored dipole-dipole matrix element values

MDQC_EXTERN int ntpair;
MDQC_EXTERN std::vector<std::vector<int>> tindex;
MDQC_EXTERN std::vector<std::vector<double>> tdipdip;
}
