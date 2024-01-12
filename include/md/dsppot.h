// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  dsppot  --  dispersion interaction scale factors  //
//                                                    //
////////////////////////////////////////////////////////

// dsp2scale   scale factor for 1-2 dispersion energy interactions
// dsp3scale   scale factor for 1-3 dispersion energy interactions
// dsp4scale   scale factor for 1-4 dispersion energy interactions
// dsp5scale   scale factor for 1-5 dispersion energy interactions
// use_dcorr   flag to use long range dispersion correction

MDQC_EXTERN real dsp2scale;
MDQC_EXTERN real dsp3scale;
MDQC_EXTERN real dsp4scale;
MDQC_EXTERN real dsp5scale;
MDQC_EXTERN bool use_dcorr;
}
