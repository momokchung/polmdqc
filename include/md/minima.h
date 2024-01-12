// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  minima  --  general parameters for minimizations  //
//                                                    //
////////////////////////////////////////////////////////

// maxiter   maximum number of iterations during optimization
// nextiter  iteration number to use for the first iteration
// fctmin    value below which function is deemed optimized
// hguess    initial value for the H-matrix diagonal elements

MDQC_EXTERN int maxiter;
MDQC_EXTERN int nextiter;
MDQC_EXTERN real fctmin;
MDQC_EXTERN real hguess;
}
