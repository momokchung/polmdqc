// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  kdsp  --  damped dispersion forcefield parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// dspsix   C6 dispersion coefficient for each atom class
// dspdmp   alpha dispersion parameter for each atom class

MDQC_EXTERN MDQCArray<real> dspsix;
MDQC_EXTERN MDQCArray<real> dspdmp;
}
