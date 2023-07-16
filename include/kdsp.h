///////////////////////////////////////////////////////////
//                                                       //
//  kdsp.h  --  damped dispersion forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// dspsix   C6 dispersion coefficient for each atom class
// dspdmp   alpha dispersion parameter for each atom class


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<double> dspsix;
QCMD_EXTERN std::vector<double> dspdmp;
