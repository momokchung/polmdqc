//////////////////////////////////////////////////////////
//                                                      //
//  scales.h  --  optimization parameter scale factors  //
//                                                      //
//////////////////////////////////////////////////////////

// scale      multiplicative factor for each optimization parameter
// set_scale  logical flag to show if scale factors have been set


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<double> scale;
QCMD_EXTERN bool set_scale;
