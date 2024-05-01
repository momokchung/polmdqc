// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "precision.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  alforder  --  reorder/partition atoms for AlphaMol2  //
//                                                       //
///////////////////////////////////////////////////////////

void alforder(real xmin, real ymin, real zmin, real xmax, real ymax, real zmax, real rmax, int nthreads, std::vector<int>& Nval);
}
