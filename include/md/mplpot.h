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

// c     PenTyp       penetration enum class (None, Gordon1, Gordon2)
// c     m2scale      scale factor for 1-2 multipole energy interactions
// c     m3scale      scale factor for 1-3 multipole energy interactions
// c     m4scale      scale factor for 1-4 multipole energy interactions
// c     m5scale      scale factor for 1-5 multipole energy interactions
// c     use_chgpen   flag to use charge penetration damped potential
// c     pentyps      string to hold penetration type
// c     pentype      type of penetration damping

enum class PenTyp
{
   None,
   Gordon1,
   Gordon2,
};

MDQC_EXTERN real m2scale,m3scale;
MDQC_EXTERN real m4scale,m5scale;
MDQC_EXTERN bool use_chgpen;
MDQC_EXTERN std::string pentyps;
MDQC_EXTERN PenTyp pentype;
}
