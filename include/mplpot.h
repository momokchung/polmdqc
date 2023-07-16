///////////////////////////////////////////////////////
//                                                   //
//  mplpot.h  --  multipole functional form details  //
//                                                   //
///////////////////////////////////////////////////////

// c     m2scale      scale factor for 1-2 multipole energy interactions
// c     m3scale      scale factor for 1-3 multipole energy interactions
// c     m4scale      scale factor for 1-4 multipole energy interactions
// c     m5scale      scale factor for 1-5 multipole energy interactions
// c     use_chgpen   flag to use charge penetration damped potential
// c     pentyp       type of penetration damping (NONE, GORDON1, GORDON2)


#pragma once
#include "macro.h"
#include <string>

QCMD_EXTERN double m2scale,m3scale;
QCMD_EXTERN double m4scale,m5scale;
QCMD_EXTERN bool use_chgpen;
QCMD_EXTERN std::string pentyp;
