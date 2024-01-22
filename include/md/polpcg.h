// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  polpcg  --  induced dipoles via the PCG solver  //
//                                                  //
//////////////////////////////////////////////////////

// mindex    index into preconditioner inverse for PCG solver
// pcgpeek   value of acceleration factor for PCG peek step
// minv      preconditioner inverse for induced dipole PCG solver
// pcgprec   flag to allow use of preconditioner with PCG solver
// pcgguess  flag to use initial PCG based on direct field

MDQC_EXTERN MDQCArray<int> mindex;
MDQC_EXTERN real pcgpeek;
MDQC_EXTERN MDQCArray<real> minv;
MDQC_EXTERN bool pcgprec;
MDQC_EXTERN bool pcgguess;
}
