// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  vdwpot  --  van der Waals functional form details  //
//                                                     //
/////////////////////////////////////////////////////////

// igauss      coefficients of Gaussian fit to vdw potential
// ngauss      number of Gaussians used in fit to vdw potential
// abuck       value of "A" constant in Buckingham vdw potential
// bbuck       value of "B" constant in Buckingham vdw potential
// cbuck       value of "C" constant in Buckingham vdw potential
// ghal        value of "gamma" in buffered 14-7 vdw potential
// dhal        value of "delta" in buffered 14-7 vdw potential
// v2scale     factor by which 1-2 vdw interactions are scaled
// v3scale     factor by which 1-3 vdw interactions are scaled
// v4scale     factor by which 1-4 vdw interactions are scaled
// v5scale     factor by which 1-5 vdw interactions are scaled
// use_vcorr   flag to use long range van der Waals correction
// vdwindex    indexing mode (atom type or class) for vdw parameters
// vdwtyp      type of van der Waals potential energy function
// radtyp      type of parameter (sigma or R-min) for atomic size
// radsiz      atomic size provided as radius or diameter
// radrule     combining rule for atomic size parameters
// epsrule     combining rule for vdw well depth parameters
// gausstyp    type of Gaussian fit to van der Waals potential

constexpr int maxgauss = 10;
MDQC_EXTERN int ngauss;
MDQC_EXTERN double igauss[maxgauss][2];
MDQC_EXTERN double abuck,bbuck,cbuck;
MDQC_EXTERN double ghal,dhal;
MDQC_EXTERN double v2scale,v3scale;
MDQC_EXTERN double v4scale,v5scale;
MDQC_EXTERN bool use_vcorr;
MDQC_EXTERN std::string vdwindex;
MDQC_EXTERN std::string radtyp;
MDQC_EXTERN std::string radsiz,gausstyp;
MDQC_EXTERN std::string radrule,epsrule;
MDQC_EXTERN std::string vdwtyp;
}
