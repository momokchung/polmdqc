///////////////////////////////////////////
//                                       //
//  calcMode  --  calculation mode type  //
//                                       //
///////////////////////////////////////////


#pragma once
#include "macro.h"

enum class calcMode
{
   NONE,
   ENERGY,
   ANALYSIS,
   VIRIAL,
   GRADIENT,
   HESSIAN,
};
