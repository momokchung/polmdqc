////////////////////////////////////////////////////////////
//                                                        //
//  extfld.h  --  applied external electric field vector  //
//                                                        //
////////////////////////////////////////////////////////////

// exfld       components of applied external electric field
// use_exfld   flag to include applied external electric field


#pragma once
#include "macro.h"

QCMD_EXTERN double exfld[3];
QCMD_EXTERN bool use_exfld;
