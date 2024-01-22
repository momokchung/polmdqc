// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "sizes.h"

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  atmlst  --  bond and angle local geometry indices  //
//                                                     //
/////////////////////////////////////////////////////////

// bndlist   numbers of the bonds involving each atom
// anglist   numbers of the angles centered on each atom
// balist    numbers of the bonds comprising each angle

MDQC_EXTERN MDQCArray2D<int,maxval> bndlist;
MDQC_EXTERN MDQCArray2D<int,maxval*(maxval-1)/2> anglist;
MDQC_EXTERN MDQCArray2D<int,2> balist;
}
