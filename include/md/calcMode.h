// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
///////////////////////////////////////////
//                                       //
//  calcMode  --  calculation mode type  //
//                                       //
///////////////////////////////////////////

enum class calcMode
{
   NONE,
   ENERGY,
   ANALYSIS,
   VIRIAL,
   GRADIENT,
   HESSIAN,
};
}
