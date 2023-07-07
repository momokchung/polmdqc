//////////////////////////////////////////////////////////
//                                                      //
//  minima.h  --  general parameters for minimizations  //
//                                                      //
//////////////////////////////////////////////////////////

// maxiter   maximum number of iterations during optimization
// nextiter  iteration number to use for the first iteration
// fctmin    value below which function is deemed optimized
// hguess    initial value for the H-matrix diagonal elements


#pragma once
#include "macro.h"

QCMD_EXTERN int maxiter;
QCMD_EXTERN int nextiter;
QCMD_EXTERN double fctmin;
QCMD_EXTERN double hguess;
