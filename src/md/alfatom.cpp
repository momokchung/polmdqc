// Author: Moses KJ Chung
// Year:   2024

#include "alfatom.h"
#include "alfp.h"
#include <cmath>

namespace polmdqc
{
/////////////////////////////////////////////
//                                         //
//  alfatom  --  atom class for AlphaMol2  //
//                                         //
/////////////////////////////////////////////

// "alfatom" class characterizes the atoms
// used in Delaunay/Alpha complex theory

AlfAtom::AlfAtom(int idx, real x, real y, real z, real r, real coefs, real coefv, real coefm, real coefg)
{
    this->index = idx;
    this->coord[0] = truncate_real(x,alfdigit);
    this->coord[1] = truncate_real(y,alfdigit);
    this->coord[2] = truncate_real(z,alfdigit);
    this->r = truncate_real(r,alfdigit);
    this->coefs = coefs;
    this->coefv = coefv;
    this->coefm = coefm;
    this->coefg = coefg;
}

AlfAtom::~AlfAtom() {}

real AlfAtom::truncate_real(real x, int ndigit)
{
    real x_out,y;
    real fact;

    int mantissa;
    int digit;

    mantissa = (int) x;
    y = x - mantissa;

    x_out = mantissa;
    fact = 1;
    for (int i = 0; i < ndigit; i++) {
        fact *= 10;
        digit = (int) std::round(y*10);
        y = 10*(y-digit/10.0);
        x_out += digit/fact;
    }

    return x_out;
}
}
