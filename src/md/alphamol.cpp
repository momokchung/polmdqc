// Author: Moses KJ Chung and Patrice Koehl
// Year:   2023

#include "alfcx.h"
#include "alfcxedges.h"
#include "alfcxfaces.h"
#include "alphamol.h"
#include "alphavol.h"
#include "alphmol.h"
#include "atomid.h"
#include "atoms.h"
#include "delaunay.h"
#include "files.h"
#include "inform.h"
#include "initdelcx.h"
#include "kvdws.h"
#include "tetrahedron.h"
#include "usage.h"
#include "vertex.h"
#include <iostream>
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  alphamol  --  compute surface area and volume  //
//                                                 //
/////////////////////////////////////////////////////

// "alphamol" computes volume, surface area, mean, and gaussian curvature
//
// literature reference:
//
// P. Koehl, A. Akopyan, and H. Edelsbrunner, "Computing the Volume,
// Surface Area, Mean, and Gaussian Curvatures of Molecules and Their
// Derivatives", Journal of Chemical Information and Modeling,
// 63, 973-985, (2023).
//
// github reference:
//
// https://github.com/pkoehl/AlphaMol

void alphamol(real r_h2o, bool deriv)
{
    clock_t start_s, stop_s;

    // perform dynamic allocation of some global arrays
    radii.allocate(n);
    coefS.allocate(n);
    coefV.allocate(n);
    coefM.allocate(n);
    coefG.allocate(n);

    // set radii and initialize coefficients
    for(int i = 0; i < n; i++) {
        if (use[i+1]) {
            radii[i] = rad[atomClass[i]] + r_h2o;
        }
        else {
            radii[i] = 0.;
        }
        coefS[i] = 1.0;
        coefV[i] = 1.0;
        coefM[i] = 1.0;
        coefG[i] = 1.0;
    }

    // initialize Delaunay procedure
    initdelcx();

    // compute Delaunay triangulation
    start_s = clock();
    delaunay();
    stop_s = clock();
    if (verbose) {
        printf("\n Delaunay compute time : %10.6e seconds\n", (stop_s-start_s)/double(CLOCKS_PER_SEC));
    }

    // generate alpha complex (with alpha=0.0)
	start_s = clock();
	real alpha = 0;
	alfcx(alpha);
	stop_s = clock();
    if (verbose) {
        printf("\n AlphaCx compute time : %10.6e seconds\n", (stop_s-start_s)/double(CLOCKS_PER_SEC));
    }

    // Compute surface area and, optionally volume of the union of balls.
    // If requested, compute also their derivatives
    alfcxedges();
    alfcxfaces();

    // allocate some dynamic arrays
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

    start_s = clock();
    alphavol(wsurf, wvol, wmean, wgauss, tsurf, tvol, tmean, tgauss,
        surf.ptr(), vol.ptr(), mean.ptr(), gauss.ptr(),
        dsurf.ptr(), dvol.ptr(), dmean.ptr(), dgauss.ptr(), deriv);
    stop_s = clock();

    if (verbose) {
        printf("\n Volumes compute time : %10.6e seconds\n", (stop_s-start_s)/double(CLOCKS_PER_SEC));
    }
}
}
