// Author: Moses KJ Chung
// Year:   2024

#include "alfatom.h"

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
    this->coord[0] = x;
    this->coord[1] = y;
    this->coord[2] = z;
    this->r = r;
    this->coefs = coefs;
    this->coefv = coefv;
    this->coefm = coefm;
    this->coefg = coefg;

    this->w = x*x + y*y + z*z - r*r;
}

AlfAtom::~AlfAtom() {}
}
