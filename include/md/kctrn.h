// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  kctrn  --  charge transfer forcefield parameters  //
//                                                    //
////////////////////////////////////////////////////////

// ctchg     charge transfer magnitude for each atom class
// ctdmp     alpha charge transfer parameter for each atom class

MDQC_EXTERN MDQCArray<real> ctchg;
MDQC_EXTERN MDQCArray<real> ctdmp;
}
