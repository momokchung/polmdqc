// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "precision.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  alfedge  --  check if edge belongs to alpha complex  //
//                                                       //
///////////////////////////////////////////////////////////

void alfedge(real* a, real* b, real ra, real rb, 
    real* cg, std::vector<int>& listcheck, int& irad, int& iattach, real alpha);
}
