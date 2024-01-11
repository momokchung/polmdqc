// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "calcMode.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  torque  --  convert single site torque to force  //
//                                                   //
///////////////////////////////////////////////////////

template <CalcMode CalculationMode>
void torque(const std::vector<std::vector<real>>* trqPtr, std::vector<std::vector<real>>* dePtr);
}
