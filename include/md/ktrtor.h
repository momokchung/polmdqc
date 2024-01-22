// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  ktrtor  --  torsion-torsion forcefield parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// maxntt    maximum number of torsion-torsion parameter entries
// maxtgrd   maximum dimension of torsion-torsion spline grid
// maxtgrd2  maximum number of torsion-torsion spline grid points
// tnx       number of columns in torsion-torsion spline grid
// tny       number of rows in torsion-torsion spline grid
// ttx       angle values for first torsion of spline grid
// tty       angle values for second torsion of spline grid
// tbf       function values at points on spline grid
// tbx       gradient over first torsion of spline grid
// tby       gradient over second torsion of spline grid
// tbxy      Hessian cross components over spline grid
// ktt       string of torsion-torsion atom classes

MDQC_EXTERN int maxntt;
constexpr int maxtgrd = 30;
constexpr int maxtgrd2= maxtgrd*maxtgrd;
MDQC_EXTERN MDQCArray<int> tnx;
MDQC_EXTERN MDQCArray<int> tny;
MDQC_EXTERN MDQCArray2D<real,maxtgrd> ttx;
MDQC_EXTERN MDQCArray2D<real,maxtgrd> tty;
MDQC_EXTERN MDQCArray2D<real,maxtgrd2> tbf;
MDQC_EXTERN MDQCArray2D<real,maxtgrd2> tbx;
MDQC_EXTERN MDQCArray2D<real,maxtgrd2> tby;
MDQC_EXTERN MDQCArray2D<real,maxtgrd2> tbxy;
MDQC_EXTERN MDQCArray<std::string> ktt;
}
