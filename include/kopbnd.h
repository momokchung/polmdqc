/////////////////////////////////////////////////////////////
//                                                         //
//  kopbnd.h  --  out-of-plane bend forcefield parameters  //
//                                                         //
/////////////////////////////////////////////////////////////

// maxnopb   maximum number of out-of-plane bending entries
// opbn      force constant parameters for out-of-plane bending
// kopb      string of atom classes for out-of-plane bending


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnopb;
QCMD_EXTERN std::vector<double> opbn;
QCMD_EXTERN std::vector<std::string> kopb;
