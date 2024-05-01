// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
///////////////////////////////////////////
//                                       //
//  alfp  --  AlphaMol2 control options  //
//                                       //
///////////////////////////////////////////

// alfmeth   algorithm to use (AlphaMol / AlphaMol2)
// alfsort   algorithm to sort and partition atoms
// alfnthd   number of threads to use for AlphaMol2

enum class AlfMethod
{
    AlphaMol,
    AlphaMol2,
};

enum class AlfSort
{
    None,
    Sort3D,
    BRIO,
    Split,
    KDTree,
};

MDQC_EXTERN AlfMethod alfmeth;
MDQC_EXTERN AlfSort alfsort;
MDQC_EXTERN int alfnthd;
}
