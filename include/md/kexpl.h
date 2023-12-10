// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  kexpl  --  exch-polarization forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// pepk     exchange-polarization spring constant for atom classes
// peppre   exchange-polarization prefactor for atom classes
// pepdmp   exchange-polarization damping alpha for atom classes
// pepl     exchange-polarization logical flag for atom classes

MDQC_EXTERN std::vector<double> pepk;
MDQC_EXTERN std::vector<double> peppre;
MDQC_EXTERN std::vector<double> pepdmp;
MDQC_EXTERN std::vector<bool> pepl;
}
