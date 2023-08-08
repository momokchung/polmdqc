///////////////////////////////////////////////////////////
//                                                       //
//  ewald.h  --  Ewald summation parameters and options  //
//                                                       //
///////////////////////////////////////////////////////////

// aewald     current value of Ewald convergence coefficient
// aeewald    Ewald convergence coefficient for electrostatics
// apewald    Ewald convergence coefficient for polarization
// adewald    Ewald convergence coefficient for dispersion
// boundary   Ewald boundary condition; none, tinfoil or vacuum


#pragma once
#include "macro.h"
#include <string>

QCMD_EXTERN double aewald;
QCMD_EXTERN double aeewald;
QCMD_EXTERN double apewald;
QCMD_EXTERN double adewald;
QCMD_EXTERN std::string boundary;
