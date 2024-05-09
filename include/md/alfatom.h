// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "precision.h"
#include <bitset>

namespace polmdqc
{
/////////////////////////////////////////////
//                                         //
//  alfatom  --  atom class for AlphaMol2  //
//                                         //
/////////////////////////////////////////////

// index   index of atom
// r       radius of vertex
// x       x Cartesian coordinate
// y       y Cartesian coordinate
// z       z Cartesian coordinate
// w       weight of vertex
// coefs   coefficient for surface
// coefv   coefficient for volume
// coefm   coefficient for mean curvature
// coefg   coefficient for gaussian curvature

class AlfAtom {
public:
    int index;
    real r;
    real coord[3];
    real coefs,coefv,coefm,coefg;

    AlfAtom() {}

    AlfAtom(int idx, real x, real y, real z, real r, real coefs, real coefv, real coefm, real coefg);

    ~AlfAtom();

private:
    real truncate_real(real x, int ndigit);
};
}
