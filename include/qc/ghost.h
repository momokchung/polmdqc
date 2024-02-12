// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  ghost  --  ghost atoms in quantum calculations  //
//                                                  //
//////////////////////////////////////////////////////

// nghst   total number of ghost atoms in the current system
// ghst    true if an atom is a ghost, false otherwise

MDQC_EXTERN int nghst;
MDQC_EXTERN MDQCArray<bool> ghst;
}
