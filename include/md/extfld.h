// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  extfld  --  applied external electric field vector  //
//                                                      //
//////////////////////////////////////////////////////////

// exfld       components of applied external electric field
// use_exfld   flag to include applied external electric field

MDQC_EXTERN real exfld[3];
MDQC_EXTERN bool use_exfld;
}
