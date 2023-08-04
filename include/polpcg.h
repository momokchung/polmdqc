////////////////////////////////////////////////////////
//                                                    //
//  polpcg.h  --  induced dipoles via the PCG solver  //
//                                                    //
////////////////////////////////////////////////////////

// mindex    index into preconditioner inverse for PCG solver
// pcgpeek   value of acceleration factor for PCG peek step
// minv      preconditioner inverse for induced dipole PCG solver
// pcgprec   flag to allow use of preconditioner with PCG solver
// pcgguess  flag to use initial PCG based on direct field


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<int> mindex;
QCMD_EXTERN double pcgpeek;
QCMD_EXTERN std::vector<double> minv;
QCMD_EXTERN bool pcgprec;
QCMD_EXTERN bool pcgguess;
