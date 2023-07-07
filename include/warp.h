//////////////////////////////////////////////////////////
//                                                      //
//  warp.h  --  potential surface smoothing parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// deform       value of the smoothing deformation parameter
// difft        diffusion coefficient for torsional potential
// diffv        diffusion coefficient for van der Waals potential
// diffc        diffusion coefficient for charge-charge potential
// m2           second moment of the GDA gaussian for each atom
// use_smooth   flag to use a potential energy smoothing method
// use_dem      flag to use diffusion equation method potential
// use_gda      flag to use gaussian density annealing potential
// use_tophat   flag to use analytical tophat smoothed potential
// use_stophat  flag to use shifted tophat smoothed potential


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN double deform;
QCMD_EXTERN double difft;
QCMD_EXTERN double diffv;
QCMD_EXTERN double diffc;
QCMD_EXTERN std::vector<double> m2;
QCMD_EXTERN bool use_smooth;
QCMD_EXTERN bool use_dem;
QCMD_EXTERN bool use_gda;
QCMD_EXTERN bool use_tophat;
QCMD_EXTERN bool use_stophat;
