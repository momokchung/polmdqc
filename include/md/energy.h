// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "calcMode.h"
#include "precision.h"
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  energy  --  evaluates energy, gradient, and virial  //
//                                                      //
//////////////////////////////////////////////////////////

template <CalcMode CalculationMode>
void energy();
}
