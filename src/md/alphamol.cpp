// Author: Moses KJ Chung and Patrice Koehl
// Year:   2023

#include "alphamol.h"
#include "alphmol.h"
#include "atomid.h"
#include "atoms.h"
#include "files.h"
#include "kvdws.h"
#include "usage.h"
#include <iostream>

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

void alphamol(real r_h2o, bool computeDeriv)
{
    // set timer
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

}
}
