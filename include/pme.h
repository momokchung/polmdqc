///////////////////////////////////////////////////////////
//                                                       //
//  pme.h  --  values for particle mesh Ewald summation  //
//                                                       //
///////////////////////////////////////////////////////////

// nfft1      current number of PME grid points along a-axis
// nfft2      current number of PME grid points along b-axis
// nfft3      current number of PME grid points along c-axis
// nefft1     number of grid points along electrostatic a-axis
// nefft2     number of grid points along electrostatic b-axis
// nefft3     number of grid points along electrostatic c-axis
// ndfft1     number of grid points along dispersion a-axis
// ndfft2     number of grid points along dispersion b-axis
// ndfft3     number of grid points along dispersion c-axis
// bsorder    current order of the PME B-spline values
// bseorder   order of the electrostatic PME B-spline values
// bsporder   order of the polarization PME B-spline values
// bsdorder   order of the dispersion PME B-spline values
// igrid      initial Ewald grid values for B-spline
// bsmod1     B-spline moduli along the a-axis direction
// bsmod2     B-spline moduli along the b-axis direction
// bsmod3     B-spline moduli along the c-axis direction
// bsbuild    B-spline derivative coefficient temporary storage
// thetai1    B-spline coefficients along the a-axis
// thetai2    B-spline coefficients along the b-axis
// thetai3    B-spline coefficients along the c-axis
// qgrid      values on the particle mesh Ewald grid
// qfac       prefactors for the particle mesh Ewald grid


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nfft1,nfft2,nfft3;
QCMD_EXTERN int nefft1,nefft2,nefft3;
QCMD_EXTERN int ndfft1,ndfft2,ndfft3;
QCMD_EXTERN int bsorder,bseorder;
QCMD_EXTERN int bsporder,bsdorder;
QCMD_EXTERN std::vector<std::vector<int>> igrid;
QCMD_EXTERN std::vector<double> bsmod1;
QCMD_EXTERN std::vector<double> bsmod2;
QCMD_EXTERN std::vector<double> bsmod3;
QCMD_EXTERN std::vector<std::vector<double>> bsbuild;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> thetai1;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> thetai2;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> thetai3;
QCMD_EXTERN std::vector<std::vector<std::vector<std::vector<double>>>> qgrid;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> qfac;
