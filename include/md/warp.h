// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  warp  --  potential surface smoothing parameters  //
//                                                    //
////////////////////////////////////////////////////////

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

MDQC_EXTERN real deform;
MDQC_EXTERN real difft;
MDQC_EXTERN real diffv;
MDQC_EXTERN real diffc;
MDQC_EXTERN std::vector<real> m2;
MDQC_EXTERN bool use_smooth;
MDQC_EXTERN bool use_dem;
MDQC_EXTERN bool use_gda;
MDQC_EXTERN bool use_tophat;
MDQC_EXTERN bool use_stophat;
}
