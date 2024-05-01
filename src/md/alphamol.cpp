// Author: Moses KJ Chung and Patrice Koehl
// Year:   2023

#include "alfcx.h"
#include "alfcxedges.h"
#include "alfcxfaces.h"
#include "alphamol.h"
#include "alphavol.h"
#include "alphmol.h"
#include "atoms.h"
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
    clock_t start_s,stop_s;
    real total = 0;

    // initialize Delaunay procedure
    start_s = clock();
    initdelcx();
    stop_s = clock();
    if (verbose) {
        printf("\n Initdelcx compute time : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    // compute Delaunay triangulation
    start_s = clock();
    delaunay();
    stop_s = clock();
    if (verbose) {
        printf("\n Delaunay compute time : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    // generate alpha complex (with alpha=0.0)
    start_s = clock();
    real alpha = 0;
    alfcx(alpha);
    stop_s = clock();
    if (verbose) {
        printf("\n AlphaCx compute time : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    // Compute surface area and, optionally volume of the union of balls.
    // If requested, compute also their derivatives
    start_s = clock();
    alfcxedges();
    alfcxfaces();
    stop_s = clock();
    if (verbose) {
        printf("\n AlphaCxEdges compute time : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    start_s = clock();
    alphavol(wsurf, wvol, wmean, wgauss,
        surf.ptr(), vol.ptr(), mean.ptr(), gauss.ptr(),
        dsurf.ptr(), dvol.ptr(), dmean.ptr(), dgauss.ptr(), deriv);
    stop_s = clock();

    if (verbose) {
        printf("\n Volumes compute time : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }
    if (verbose) {
        printf("\n AlphaMol compute time : %10.6f ms\n", total*1000);
    }
}
}
