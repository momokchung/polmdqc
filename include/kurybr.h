/////////////////////////////////////////////////////////////
//                                                         //
//  kurybr.h  --  Urey-Bradley term forcefield parameters  //
//                                                         //
/////////////////////////////////////////////////////////////

// maxnu   maximum number of Urey-Bradley parameter entries
// ucon    force constant parameters for Urey-Bradley terms
// dst13   ideal 1-3 distance parameters for Urey-Bradley terms
// ku      string of atom classes for Urey-Bradley terms


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnu;
QCMD_EXTERN std::vector<double> ucon;
QCMD_EXTERN std::vector<double> dst13;
QCMD_EXTERN std::vector<std::string> ku;
