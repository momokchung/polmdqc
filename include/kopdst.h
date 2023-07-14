/////////////////////////////////////////////////////////////
//                                                         //
//  kopdst.h  --  out-of-plane distance forcefield params  //
//                                                         //
/////////////////////////////////////////////////////////////

// maxnopd   maximum number of out-of-plane distance entries
// opds      force constant parameters for out-of-plane distance
// kopd      string of atom classes for out-of-plane distance


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnopd;
QCMD_EXTERN std::vector<double> opds;
QCMD_EXTERN std::vector<std::string> kopd;
