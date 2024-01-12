// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  linmin  --  line search minimization parameters  //
//                                                   //
///////////////////////////////////////////////////////

// intmax   maximum number of interpolations during line search
// stpmin   minimum step length in current line search direction
// stpmax   maximum step length in current line search direction
// cappa    stringency of line search (0=tight < cappa < 1=loose)
// slpmax   projected gradient above which stepsize is reduced
// angmax   maximum angle between search direction and -gradient

MDQC_EXTERN int intmax;
MDQC_EXTERN real stpmin;
MDQC_EXTERN real stpmax;
MDQC_EXTERN real cappa;
MDQC_EXTERN real slpmax;
MDQC_EXTERN real angmax;
}
