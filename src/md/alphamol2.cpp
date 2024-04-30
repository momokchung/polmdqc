// Author: Moses KJ Chung
// Year:   2024

#include "alphamol2.h"
#include "alphmol.h"
#include "atoms.h"
#include "dlauny2.h"
#include "hilbert.h"
#include "initalfatm.h"
#include "openmp.h"
#include "retmax.h"
#include "retmin.h"
#include <algorithm>
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  alphamol2  --  parallel surface area and volume  //
//                                                   //
///////////////////////////////////////////////////////

// "alphamol2" computes volume, surface area, mean, and
// gaussian curvature in parallel
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

void alphamol2(bool deriv)
{
    clock_t start_s, stop_s;

    // initialize the alfatoms vector
    initalfatm();

    // if needed, reorder  atoms
    real xmin,ymin,zmin;
    real xmax,ymax,zmax;
    real rmax;

    xmin = retmin(x);
    ymin = retmin(y);
    zmin = retmin(z);
    xmax = retmax(x);
    ymax = retmax(y);
    zmax = retmax(z);
    rmax = retmax(radii);

    real t1,t2,u1,u2,diff;
    std::vector<int> Nval(nthread + 1, 0);

    printf("nthread %d\n", nthread);
}
}
