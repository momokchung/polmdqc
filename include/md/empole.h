// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "calcMode.h"
#include "mplpot.h"

namespace polmdqc
{
////////////////////////////////////////////////
//                                            //
//  empole  --  atomic multipole calculation  //
//                                            //
////////////////////////////////////////////////

template <CalcMode CalculationMode>
void empole();

template <CalcMode CalculationMode, PenTyp PenType>
void empole_a();

// template <CalcMode CalculationMode>
// void empole_b();

// template <CalcMode CalculationMode>
// void empole_c();

// template <CalcMode CalculationMode>
// void empole_d();
}
