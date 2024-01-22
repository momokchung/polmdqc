// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

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

MDQC_EXTERN MDQCArray<real> pepk;
MDQC_EXTERN MDQCArray<real> peppre;
MDQC_EXTERN MDQCArray<real> pepdmp;
MDQC_EXTERN MDQCArray<bool> pepl;
}
