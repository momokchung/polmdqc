// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  reppot  --  repulsion functional form details  //
//                                                 //
/////////////////////////////////////////////////////

// r2scale    scale factor for 1-2 repulsion energy interactions
// r3scale    scale factor for 1-3 repulsion energy interactions
// r4scale    scale factor for 1-4 repulsion energy interactions
// r5scale    scale factor for 1-5 repulsion energy interactions

MDQC_EXTERN real r2scale;
MDQC_EXTERN real r3scale;
MDQC_EXTERN real r4scale;
MDQC_EXTERN real r5scale;
}
