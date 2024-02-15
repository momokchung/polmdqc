// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////////////
//                                                          //
//  cheb  --  store parameters for Chebyshev interpolation  //
//                                                          //
//////////////////////////////////////////////////////////////

// intorder    order of Chebyshev interpolation
// chebmmax    maximum m value for evaluating Boys function
// chebnint    number of Chebyshev intervals
// mmax        maximum order to compute reference
// chebtmax    maximum t value for evaluating Boys function
// chebdlta    width of Chebyshev interval
// chebdlta2   half width of Chebyshev interval
// ihalf       store half values (-0.5, 0.5, 1.5, 2.5, ...)
// chebtble    table for Chebyshev coefficients

constexpr int intorder = 7;
constexpr int chebmmax = 100;
constexpr int chebnint = 819;
constexpr int mmax = 200;
constexpr double chebtmax = 117;
constexpr double chebdlta = 1./7;
constexpr double chebdlta2 = chebdlta/2;
MDQC_EXTERN double ihalf[mmax];
MDQC_EXTERN double chebtble[chebnint][(chebmmax+1)*(intorder+1)];
}
