// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  chgpot  --  charge-charge functional form details  //
//                                                     //
/////////////////////////////////////////////////////////

// electric   energy factor in kcal/mole for current force field
// dielec     dielectric constant for electrostatic interactions
// ebuffer    electrostatic buffering constant added to distance
// c1scale    factor by which 1-1 charge interactions are scaled
// c2scale    factor by which 1-2 charge interactions are scaled
// c3scale    factor by which 1-3 charge interactions are scaled
// c4scale    factor by which 1-4 charge interactions are scaled
// c5scale    factor by which 1-5 charge interactions are scaled
// neutnbr    logical flag governing use of neutral group neighbors
// neutcut    logical flag governing use of neutral group cutoffs

MDQC_EXTERN real electric;
MDQC_EXTERN real dielec,ebuffer;
MDQC_EXTERN real c1scale,c2scale;
MDQC_EXTERN real c3scale,c4scale;
MDQC_EXTERN real c5scale;
MDQC_EXTERN bool neutnbr,neutcut;
}
