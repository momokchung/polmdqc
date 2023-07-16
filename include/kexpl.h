////////////////////////////////////////////////////////////
//                                                        //
//  kexpl.h  --  exch-polarization forcefield parameters  //
//                                                        //
////////////////////////////////////////////////////////////

// pepk     exchange-polarization spring constant for atom classes
// peppre   exchange-polarization prefactor for atom classes
// pepdmp   exchange-polarization damping alpha for atom classes
// pepl     exchange-polarization logical flag for atom classes


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<double> pepk;
QCMD_EXTERN std::vector<double> peppre;
QCMD_EXTERN std::vector<double> pepdmp;
QCMD_EXTERN std::vector<bool> pepl;
