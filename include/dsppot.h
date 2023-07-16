//////////////////////////////////////////////////////////
//                                                      //
//  dsppot.h  --  dispersion interaction scale factors  //
//                                                      //
//////////////////////////////////////////////////////////

// dsp2scale   scale factor for 1-2 dispersion energy interactions
// dsp3scale   scale factor for 1-3 dispersion energy interactions
// dsp4scale   scale factor for 1-4 dispersion energy interactions
// dsp5scale   scale factor for 1-5 dispersion energy interactions
// use_dcorr   flag to use long range dispersion correction


#pragma once
#include "macro.h"

QCMD_EXTERN double dsp2scale;
QCMD_EXTERN double dsp3scale;
QCMD_EXTERN double dsp4scale;
QCMD_EXTERN double dsp5scale;
QCMD_EXTERN bool use_dcorr;
