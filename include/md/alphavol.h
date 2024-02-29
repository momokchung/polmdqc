// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "dlauny.h"
#include "precision.h"

namespace polmdqc
{
///////////////////////////////////////////////////////////////
//                                                           //
//  alphavol  --  compute surf, vol, mean and gaussian curv  //
//                                                           //
///////////////////////////////////////////////////////////////

void alphavol(real& WSurf, real& WVol, real& WMean, real& WGauss,
    real& Surf, real& Vol, real& Mean, real& Gauss,
    real* ballwsurf, real* ballwvol, real* ballwmean, real* ballwgauss,
    real* dsurf_coord, real* dvol_coord, real* dmean_coord, real* dgauss_coord, bool compder);
}
