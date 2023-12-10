// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  mplpot  --  multipole functional form details  //
//                                                 //
/////////////////////////////////////////////////////

// c     m2scale      scale factor for 1-2 multipole energy interactions
// c     m3scale      scale factor for 1-3 multipole energy interactions
// c     m4scale      scale factor for 1-4 multipole energy interactions
// c     m5scale      scale factor for 1-5 multipole energy interactions
// c     use_chgpen   flag to use charge penetration damped potential
// c     pentyp       type of penetration damping (NONE, GORDON1, GORDON2)

MDQC_EXTERN double m2scale,m3scale;
MDQC_EXTERN double m4scale,m5scale;
MDQC_EXTERN bool use_chgpen;
MDQC_EXTERN std::string pentyp;
}
