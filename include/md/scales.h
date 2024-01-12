// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  scales  --  optimization parameter scale factors  //
//                                                    //
////////////////////////////////////////////////////////

// scale      multiplicative factor for each optimization parameter
// set_scale  logical flag to show if scale factors have been set

MDQC_EXTERN std::vector<real> scale;
MDQC_EXTERN bool set_scale;
}
