//////////////////////////////////////////////////////////
//                                                      //
//  kctrn.h  --  charge transfer forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// ctchg     charge transfer magnitude for each atom class
// ctdmp     alpha charge transfer parameter for each atom class


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<double> ctchg;
QCMD_EXTERN std::vector<double> ctdmp;
