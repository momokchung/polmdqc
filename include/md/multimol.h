// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "precision.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  multimol  --  launch multithreaded AlphaMol  //
//                                               //
///////////////////////////////////////////////////

void multimol(real buffer, bool deriv, int nthreads, std::vector<int>& Nval);
}
