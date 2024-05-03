// Author: Moses KJ Chung
// Year:   2024

#include "alfatom.h"
#include "alphmol.h"
#include "atoms.h"
#include "dlauny2.h"
#include <algorithm>

namespace polmdqc
{
//////////////////////////////////////////////////
//                                              //
//  initalfatm  --  initialize AlphaMol2 Atoms  //
//                                              //
//////////////////////////////////////////////////

// "initalfatm" initializes the alfatoms vector

void initalfatm()
{
    // resize alfatoms
    alfatoms.clear();
    alfatoms.reserve(n);

    // copy atoms into alfatoms list
    real xi, yi, zi, ri, cs, cv, cm, cg;
    for(int i = 0; i < n; i++) {
        xi = x[i];
        yi = y[i];
        zi = z[i];
        ri = radii[i];
        cs = coefS[i];
        cv = coefV[i];
        cm = coefM[i];
        cg = coefG[i];
        AlfAtom atm(i, xi, yi, zi, ri, cs, cv, cm, cg);
        alfatoms.push_back(atm);
    }

    // If needed, randomly shuffle the atoms
    bool shuffle = false;
    if (shuffle) {
        std::random_shuffle(alfatoms.begin(), alfatoms.end());
    }
}
}
