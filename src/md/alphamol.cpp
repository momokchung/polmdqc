// Author: Moses KJ Chung and Patrice Koehl
// Year:   2023

#include "alfcx.h"
#include "alfcxedges.h"
#include "alfcxfaces.h"
#include "alphamol.h"
#include "alphavol.h"
#include "alphmol.h"
#include "delaunay.h"
#include "files.h"
#include "inform.h"
#include "initdelcx.h"
#include "tetrahedron.h"
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

// "alphamol" computes volume, surface area, mean, and
// gaussian curvature
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

void alphamol(bool deriv)
{
    clock_t start_s, stop_s;

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
