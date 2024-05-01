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

// wsurf     weighted surface area of the system
// wvol      weighted volume of the system
// wmean     weighted mean curvature of the system
// wgauss    weighted gaussian curvature of the system
// coord     coordinates for each atom
// radii     radius for each atom
// coefS     coefficient for surface area calculation
// coefV     coefficient for volume calculation
// coefM     coefficient for mean curvature calculation
// coefG     coefficient for gaussian curvature calculation
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

MDQC_EXTERN real wsurf;
MDQC_EXTERN real wvol;
MDQC_EXTERN real wmean;
MDQC_EXTERN real wgauss;
MDQC_EXTERN MDQCArray<real> radii;
MDQC_EXTERN MDQCArray<real> coefS;
MDQC_EXTERN MDQCArray<real> coefV;
MDQC_EXTERN MDQCArray<real> coefM;
MDQC_EXTERN MDQCArray<real> coefG;
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
