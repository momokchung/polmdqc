// Author: Moses KJ Chung
// Year:   2023

#pragma once

#include "precision.h"
#include <cmath>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  libfunc  --  math functions for float or double  //
//                                                   //
///////////////////////////////////////////////////////

#if POLMDQC_REAL_SIZE == 8
#   define REAL_SQRT     sqrt
#   define REAL_EXP      exp
#   define REAL_FLOOR    floor
#   define REAL_POW      pow
#   define REAL_COS      cos
#   define REAL_SIN      sin
#   define REAL_ACOS     acos
#   define REAL_ASIN     asin
#   define REAL_SINH     sinh
#   define REAL_COSH     cosh
#   define REAL_ERF      erf
#   define REAL_ERFC(x)  (1 - erf(x))
#   define REAL_MIN      fmin
#   define REAL_MAX      fmax
#   define REAL_SIGN     copysign
#   define REAL_ABS      std::fabs
#endif
}
