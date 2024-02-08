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

#if POLMDQC_REALQ_SIZE == 8
#   define REALQ_SQRT     sqrt
#   define REALQ_EXP      exp
#   define REALQ_FLOOR    floor
#   define REALQ_POW      pow
#   define REALQ_COS      cos
#   define REALQ_SIN      sin
#   define REALQ_ACOS     acos
#   define REALQ_ASIN     asin
#   define REALQ_SINH     sinh
#   define REALQ_COSH     cosh
#   define REALQ_ERF      erf
#   define REALQ_ERFC(x)  (1 - erf(x))
#   define REALQ_MIN      fmin
#   define REALQ_MAX      fmax
#   define REALQ_SIGN     copysign
#   define REALQ_ABS      std::fabs
#endif
}
