/////////////////////////////////////////////////////////
//                                                     //
//  kpolpr.h  --  special Thole forcefield parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// maxnpp   maximum number of special pair polarization entries
// thlpr    Thole damping values for special polarization pairs
// thdpr    Thole direct damping for special polarization pairs
// kppr     string of atom types for special polarization pairs


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnpp;
QCMD_EXTERN std::vector<double> thlpr;
QCMD_EXTERN std::vector<double> thdpr;
QCMD_EXTERN std::vector<std::string> kppr;
