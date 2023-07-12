///////////////////////////////////////////////////////////
//                                                       //
//  angbnd.h  --  bond angle bends in current structure  //
//                                                       //
///////////////////////////////////////////////////////////

// nangle   total number of angle bends in the system
// iang     numbers of the atoms in each angle bend
// ak       harmonic angle force constant (kcal/mole/rad**2)
// anat     ideal bond angle or phase shift angle (degrees)
// afld     periodicity for Fourier angle bending term


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nangle;
QCMD_EXTERN std::vector<std::vector<int>> iang;
QCMD_EXTERN std::vector<double> ak;
QCMD_EXTERN std::vector<double> anat;
QCMD_EXTERN std::vector<double> afld;
