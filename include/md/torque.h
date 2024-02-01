// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "calcMode.h"

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  torque  --  convert single site torque to force  //
//                                                   //
///////////////////////////////////////////////////////

template <CalcMode CalculationMode>
void torque(const real* trq, real* de, real (&virial)[3][3]);
}
