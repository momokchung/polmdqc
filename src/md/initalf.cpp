// Author: Moses KJ Chung
// Year:   2024

#include "alphmol.h"
#include "atomid.h"
#include "atoms.h"
#include "initalf.h"
#include "kvdws.h"
#include "usage.h"

namespace polmdqc
{
////////////////////////////////////////
//                                    //
//  initalf  --  initialize AlphaMol  //
//                                    //
////////////////////////////////////////

// "initalf" initialize variables used in AlphaMol

void initalf(real scoef, real vcoef, real exclude, bool deriv)
{
    // perform dynamic allocation of some global arrays
    radii.allocate(n);
    coefS.allocate(n);
    coefV.allocate(n);
    coefM.allocate(n);
    coefG.allocate(n);

    // perform dynamic allocation of some global arrays
    int nfudge = 8;
    surf.allocate(n+nfudge);
    vol.allocate(n+nfudge);
    mean.allocate(n+nfudge);
    gauss.allocate(n+nfudge);
    if (deriv) {
        dsurf.allocate(3*(n+nfudge));
        dvol.allocate(3*(n+nfudge));
        dmean.allocate(3*(n+nfudge));
        dgauss.allocate(3*(n+nfudge));
    }

    // set radii and initialize coefficients
    for(int i = 0; i < n; i++) {
        if (use[i+1]) {
            radii[i] = rad[atomClass[i]] + exclude;
        }
        else {
            radii[i] = 0.;
        }
        coefS[i] = scoef;
        coefV[i] = vcoef;
        coefM[i] = 1.0;
        coefG[i] = 1.0;
    }
}
}
