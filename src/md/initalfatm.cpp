// Author: Moses KJ Chung
// Year:   2024

#include "alfatom.h"
#include "alphmol.h"
#include "atomid.h"
#include "atoms.h"
#include "chksymm.h"
#include "dlauny2.h"
#include "random.h"
#include <algorithm>

namespace polmdqc
{
/////////////////////////////////////////////////
//                                             //
//  initalfatm  --  initialize AlphaMol Atoms  //
//                                             //
/////////////////////////////////////////////////

// "initalfatm" initializes the alfatoms vector

void initalfatm(bool deriv)
{
    // resize alfatoms
    alfatoms.clear();
    alfatoms.reserve(n);

    // initialize atomic
    for (int i = 0; i < n; i++) {
        surf[i] = 0;
        vol[i] = 0;
        mean[i] = 0;
        gauss[i] = 0;
    }

    // initialize derivatives
    if (deriv) {
        for (int i = 0; i < 3*n; i++) {
            dsurf[i] = 0;
            dvol[i] = 0;
            dmean[i] = 0;
            dgauss[i] = 0;
        }
    }

    // check for symmetry
    SymTyp symtyp;

    // allocate
    real* xref = new real[n];
    real* yref = new real[n];
    real* zref = new real[n];

    // copy coordinates
    for (int i = 0; i < n; i++) {
        xref[i] = x[i];
        yref[i] = y[i];
        zref[i] = z[i];
    }
    chksymm(n, mass.ptr(), xref, yref, zref, symtyp);

    // wiggle if system is symmetric
    real eps = 1e-5;
    bool dowiggle_linear = symtyp == SymTyp::Linear and n > 2;
    bool dowiggle_planar = symtyp == SymTyp::Planar and n > 3;
    bool dowiggle = dowiggle_linear or dowiggle_planar;
    if (dowiggle) wiggle(n, xref, yref, zref, eps);

    // copy atoms into alfatoms list
    real xi, yi, zi, ri, cs, cv, cm, cg;
    for (int i = 0; i < n; i++) {
        ri = radii[i];
        if (ri == 0) continue;
        xi = xref[i];
        yi = yref[i];
        zi = zref[i];
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

    // deallocate
    delete[] xref;
    delete[] yref;
    delete[] zref;
}
}
