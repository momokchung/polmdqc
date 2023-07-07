/////////////////////////////////////////////////////////
//                                                     //
//  linmin.h  --  line search minimization parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// intmax   maximum number of interpolations during line search
// stpmin   minimum step length in current line search direction
// stpmax   maximum step length in current line search direction
// cappa    stringency of line search (0=tight < cappa < 1=loose)
// slpmax   projected gradient above which stepsize is reduced
// angmax   maximum angle between search direction and -gradient


#pragma once
#include "macro.h"

QCMD_EXTERN int intmax;
QCMD_EXTERN double stpmin;
QCMD_EXTERN double stpmax;
QCMD_EXTERN double cappa;
QCMD_EXTERN double slpmax;
QCMD_EXTERN double angmax;
