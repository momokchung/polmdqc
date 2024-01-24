// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  alphmol  --  surface, volume, and curvature  //
//                                               //
///////////////////////////////////////////////////

// tsurf     total surface area of the system
// tvol      total volume of the system
// tmean     total mean curvature of the system
// tgauss    total gaussian curvature of the system
// wsurf     weighted surface area of the system
// wvol      weighted volume of the system
// wmean     weighted mean curvature of the system
// wgauss    weighted gaussian curvature of the system
// surf      surface area for each atom
// vol       volume for each atom
// mean      mean curvature for each atom
// gauss     gaussian curvature for each atom
// dsurf     surface area analytical gradient
// dvol      volume analytical gradient
// dmean     mean curvature analytical gradient
// dgauss    gaussian curvature analytical gradient
// ndsurf    surface area cartesian numerical gradient
// ndvol     volume cartesian numerical gradient
// ndmean    mean curvature cartesian numerical gradient
// ndgauss   gaussian curvature cartesian numerical gradient

MDQC_EXTERN real tsurf;
MDQC_EXTERN real tvol;
MDQC_EXTERN real tmean;
MDQC_EXTERN real tgauss;
MDQC_EXTERN real wsurf;
MDQC_EXTERN real wvol;
MDQC_EXTERN real wmean;
MDQC_EXTERN real wgauss;
MDQC_EXTERN MDQCArray<real> surf;
MDQC_EXTERN MDQCArray<real> vol;
MDQC_EXTERN MDQCArray<real> mean;
MDQC_EXTERN MDQCArray<real> gauss;
MDQC_EXTERN MDQCArray<real> dsurf;
MDQC_EXTERN MDQCArray<real> dvol;
MDQC_EXTERN MDQCArray<real> dmean;
MDQC_EXTERN MDQCArray<real> dgauss;
MDQC_EXTERN MDQCArray<real> ndsurf;
MDQC_EXTERN MDQCArray<real> ndvol;
MDQC_EXTERN MDQCArray<real> ndmean;
MDQC_EXTERN MDQCArray<real> ndgauss;
}
