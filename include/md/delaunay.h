// Author: Moses KJ Chung
// Year:   2024

#pragma once

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  delaunay  --  weighted Delaunay triangulation  //
//                                                 //
/////////////////////////////////////////////////////

void delaunay();

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
