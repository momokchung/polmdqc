// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "precision.h"

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  torphase  --  torsional amplitude and phase  //
//                                               //
///////////////////////////////////////////////////

void torphase(int (&ft)[6], real (&vt)[6], real (&st)[6]);
}
