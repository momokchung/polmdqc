// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "calcMode.h"

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  dcflux  --  charge flux gradient chain rule  //
//                                               //
///////////////////////////////////////////////////

template <CalcMode CalculationMode>
void dcflux(const real* cfpot, real* de, real (&virial)[3][3]);
}
