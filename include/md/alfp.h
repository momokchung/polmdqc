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

// alfmeth    algorithm to use (AlphaMol / AlphaMol2)
// alfsort    algorithm to sort and partition atoms
// alfdigit   number of digits to store for AlphaMol x,y,z,r
// alfnthd    number of threads to use for AlphaMol2
// delcxeps   sos gmp eps value for delcx
// alfcxeps   sos gmp eps value for alfcx
// alfsos     flag to use multiprecision simulation of simplicity

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
MDQC_EXTERN int alfdigit;
MDQC_EXTERN int alfnthd;
MDQC_EXTERN real delcxeps;
MDQC_EXTERN real alfcxeps;
MDQC_EXTERN bool alfsos;
}
