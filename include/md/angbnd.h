// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

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
MDQC_EXTERN std::vector<std::vector<int>> iang;
MDQC_EXTERN std::vector<double> ak;
MDQC_EXTERN std::vector<double> anat;
MDQC_EXTERN std::vector<double> afld;
}
