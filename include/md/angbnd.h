// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  angbnd  --  bond angle bends in current structure  //
//                                                     //
/////////////////////////////////////////////////////////

// nangle   total number of angle bends in the system
// iang     numbers of the atoms in each angle bend
// ak       harmonic angle force constant (kcal/mole/rad**2)
// anat     ideal bond angle or phase shift angle (degrees)
// afld     periodicity for Fourier angle bending term

MDQC_EXTERN int nangle;
MDQC_EXTERN MDQCArray2D<int,4> iang;
MDQC_EXTERN MDQCArray<real> ak;
MDQC_EXTERN MDQCArray<real> anat;
MDQC_EXTERN MDQCArray<real> afld;
}
