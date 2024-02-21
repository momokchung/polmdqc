/* ===============================================================================================
   AlphaMol: a program for computing geometric measures of a union of balls

   Author:  Patrice Koehl  (collaboration with Herbert Edelsbrunner
   Date:    9/22/2019
   Version: 1
   =============================================================================================== */

#pragma once
#include "precision.h"

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  alphamol  --  compute surface area and volume  //
//                                                 //
/////////////////////////////////////////////////////

void alphamol(real r_h2o, bool computeDeriv);
}
