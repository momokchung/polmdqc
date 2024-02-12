// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  scft  --  scf procedure settings for quantum  //
//                                                //
////////////////////////////////////////////////////

// SCFType   enum class for scf algorithm for electron repulsion integrals
// scftyp    scf algorithm for electron repulsion integrals


enum class SCFType
{
    PK,
    DF,
    Direct,
};
MDQC_EXTERN SCFType scftyp;
}
