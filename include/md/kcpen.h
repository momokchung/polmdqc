// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  kcpen  --  charge penetration forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// cpele     valence electron magnitude for each atom class
// cpalp     alpha charge penetration parameter for each atom class

MDQC_EXTERN MDQCArray<real> cpele;
MDQC_EXTERN MDQCArray<real> cpalp;
}
