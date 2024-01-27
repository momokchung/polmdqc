// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  cflux  --  charge flux terms in current structure  //
//                                                     //
/////////////////////////////////////////////////////////

// nbflx    total number of bond charge flux interactions
// naflx    total number of angle charge flux interactions
// bflx     bond stretching charge flux constant (electrons/Ang)
// aflx     angle bending charge flux constant (electrons/radian)
// abflx    asymmetric stretch charge flux constant (electrons/Ang)
// pdelta   change in partial charge
// pot      potential array for charge flux

MDQC_EXTERN int nbflx;
MDQC_EXTERN int naflx;
MDQC_EXTERN MDQCArray<real> bflx;
MDQC_EXTERN MDQCArray<real> pdelta;
MDQC_EXTERN MDQCArray<real> pot;
MDQC_EXTERN MDQCArray2D<real,2> aflx;
MDQC_EXTERN MDQCArray2D<real,2> abflx;
}
