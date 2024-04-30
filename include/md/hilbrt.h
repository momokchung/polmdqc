// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
///////////////////////////////////////////
//                                       //
//  hilbrt  --  Hilbert curve variables  //
//                                       //
///////////////////////////////////////////

constexpr int hilbert_order = 52;
constexpr int hilbert_limit = 8;
constexpr int brio_threshold = 64;
constexpr real brio_ratio = 0.125;
MDQC_EXTERN int transgc[8][3][8];
MDQC_EXTERN int tsb1mod3[8];
}
