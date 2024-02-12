// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////
//                                            //
//  sym  --  symmetry parameters for quantum  //
//                                            //
////////////////////////////////////////////////

// Symmetry   enum class for symmetry type
// sym        symmetry for current structure

enum class Symmetry
{
    C1,
    Ci,
    C2,
    Cs,
    D2,
    C2v,
    C2h,
    D2h,
};
MDQC_EXTERN Symmetry sym;
}
