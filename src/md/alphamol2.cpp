// Author: Moses KJ Chung
// Year:   2024

#include "alforder.h"
#include "alfp.h"
#include "alphamol2.h"
#include "alphmol.h"
#include "atoms.h"
#include "dlauny2.h"
#include "hilbert.h"
#include "inform.h"
#include "initalfatm.h"
#include "retmax.h"
#include "retmin.h"
#include <algorithm>
#include <sys/time.h>
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

inline void gettime(real& t1, real& u1)
{
    timeval tim;
    gettimeofday(&tim,NULL);
    t1 = tim.tv_sec;
    u1 = tim.tv_usec;
}

inline real gettimediff(real t1, real u1, real t2, real u2)
{
    return (t2-t1) + (u2-u1)*1.e-6;;
}

void alphamol2(bool deriv)
{
    real t1,t2,u1,u2,diff;
    real total = 0;

    // initialize the alfatoms vector
    gettime(t1, u1);
    initalfatm();
    gettime(t2, u2);
    diff = gettimediff(t1, u1, t2, u2);
    if (verbose) {
        printf("\n Initalfatm compute time : %10.6f ms\n", diff*1000);
        total += diff;
    }

    // if needed, reorder  atoms
    gettime(t1, u1);
    real xmin,ymin,zmin;
    real xmax,ymax,zmax;
    real rmax;
    std::vector<int> Nval(alfnthd + 1, 0);
    xmin = retmin(x);
    ymin = retmin(y);
    zmin = retmin(z);
    xmax = retmax(x);
    ymax = retmax(y);
    zmax = retmax(z);
    rmax = retmax(radii);
    alforder(xmin,ymin,zmin,xmax,ymax,zmax,rmax,alfnthd,Nval);
    gettime(t2, u2);
    diff = gettimediff(t1, u1, t2, u2);
    if (verbose) {
        printf("\n Alforder compute time : %10.6f ms\n", diff*1000);
        total += diff;
    }

    // run AlphaMol algorithm
    gettime(t1, u1);
    int natoms = alfatoms.size();
    // double WGeom[4];
    // double *geom = new double[4*(natoms+8)];
    // double *geom_order = new double[4*(natoms+8)];
    // double *dgeom = new double[12*(natoms+8)];
    double buffer = 2*rmax;

    gettime(t2, u2);
    diff = gettimediff(t1, u1, t2, u2);
    if (verbose) {
        printf("\n MeasureMol compute time : %10.6f ms\n", diff*1000);
        total += diff;
    }

    printf("alfnthd %d\n", alfnthd);
    if (verbose) {
        printf("\n AlphaMol2 compute time : %10.6f ms\n", total*1000);
    }
}
}
