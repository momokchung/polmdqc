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
void torque(const MDQCArray2D<real,3>& trq, MDQCArray2D<real,3>& de);
}
