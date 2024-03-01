// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "libfunc.h"

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  alffunc  --  various functions used in AlphaMol  //
//                                                   //
///////////////////////////////////////////////////////

// "alffunc" contains various functions used in AlphaMol

inline real psub(real r1, real r2)
{
    constexpr real epsabs = 1e-14;
    constexpr real epsrel = 1e-14;
    real r = r1 - r2;
    if (REAL_ABS(r) < epsabs) return 0.;
    real maxr = REAL_MAX(REAL_ABS(r1),REAL_ABS(r2));
    if (maxr != 0.) {
        if (REAL_ABS(r/maxr) < epsrel) return 0.;
    }
    return r;
}

inline real padd(real r1, real r2)
{
    constexpr real epsabs = 1e-14;
    constexpr real epsrel = 1e-14;
    real r = r1 + r2;
    if (REAL_ABS(r) < epsabs) return 0.;
    real maxr = REAL_MAX(REAL_ABS(r1),REAL_ABS(r2));
    if (maxr != 0.) {
        if (REAL_ABS(r/maxr) < epsrel) return 0.;
    }
    return r;
}

// builds the weight of a point: w = x**2 + y**2 + z**2 - ra**2
inline void buildweight(real ax, real ay, real az, real r, real& w)
{
    real temp1,temp2;
    temp1 = r * r;
    temp2 = ax * ax;
    temp1 = psub(temp2,temp1);
    temp2 = ay * ay;
    temp1 = padd(temp2,temp1);
    temp2 = az * az;

    w = padd(temp2,temp1);
}

inline int sgn(real d)
{
    if (d == 0.) return 0;
    else if (d > 0.) return 1;
    else return -1;
}

// sign associated with the missing infinite point
inline void missinf_sign(int i, int j, int k, int& l, int& sign)
{
    int a, b, c, d;
    l = 6 - i - j - k;
    a = i;
    b = j;
    c = k;
    sign = 1;
    if (a > b) {
        d = a;
        a = b;
        b = d;
        sign = -sign;
    }
    if (a > c) {
        d = a;
        a = c;
        c = d;
        sign = -sign;
    }
    if (b > c) {
        sign = -sign;
    }
}
}
