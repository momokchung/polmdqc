// Author: Moses KJ Chung
// Year:   2024

#include "chksymm.h"
#include "inertia.h"
#include "libfunc.h"

namespace polmdqc
{
////////////////////////////////////////////////////////////////
//                                                            //
//  subroutine chksymm  --  test for 1D, 2D & other symmetry  //
//                                                            //
////////////////////////////////////////////////////////////////

// "chksymm" examines the current coordinates for linearity,
// planarity, an internal mirror plane or center of inversion

void chksymm(int n, real* mass, real* xref, real* yref, real* zref, SymTyp& symtyp)
{
    // copy coordinates
    real* x = new real[n];
    real* y = new real[n];
    real* z = new real[n];
    for (int i = 0; i < n; i++) {
        x[i] = xref[i];
        y[i] = yref[i];
        z[i] = zref[i];
    }

    // move the atomic coordinates into the inertial frame
    inertia(2, n, mass, x, y, z);

    real eps = 1e-4;
    symtyp = SymTyp::None;
    bool xnul = true;
    bool ynul = true;
    bool znul = true;
    for (int i = 0; i < n; i++) {
        if (REAL_ABS(x[i]) > eps) xnul = false;
        if (REAL_ABS(y[i]) > eps) ynul = false;
        if (REAL_ABS(z[i]) > eps) znul = false;
    }
    if (n == 3) symtyp = SymTyp::Planar;
    if (xnul) symtyp = SymTyp::Planar;
    if (ynul) symtyp = SymTyp::Planar;
    if (znul) symtyp = SymTyp::Planar;
    if (n == 2) symtyp = SymTyp::Linear;
    if (xnul and ynul) symtyp = SymTyp::Linear;
    if (xnul and znul) symtyp = SymTyp::Linear;
    if (ynul and znul) symtyp = SymTyp::Linear;
    if (n == 1) symtyp = SymTyp::Single;

    // test mean coords for mirror plane and inversion center
    if (symtyp == SymTyp::None) {
        real xave = 0;
        real yave = 0;
        real zave = 0;
        for (int i = 0; i < n; i++) {
            xave += x[i];
            yave += y[i];
            zave += z[i];
        }
        xave = REAL_ABS(xave) / n;
        yave = REAL_ABS(yave) / n;
        zave = REAL_ABS(zave) / n;
        int nave = 0;
        if (xave < eps) nave++;
        if (yave < eps) nave++;
        if (zave < eps) nave++;
        if (nave != 0) symtyp = SymTyp::Mirror;
        if (nave == 3) symtyp = SymTyp::Center;
    }

    // deallocate
    delete[] x;
    delete[] y;
    delete[] z;
}
}
