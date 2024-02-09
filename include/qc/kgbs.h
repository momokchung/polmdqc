// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"
#include "ptable.h"
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////
//                                           //
//  kgbs  --  gaussian basis set parameters  //
//                                           //
///////////////////////////////////////////////

// AngMom      enum class for angular momentum
// BasisType   enum class for cartesian or spherical gaussian basis type
// maxnprim    maximum number of primitives for a gaussian basis function
// nprimgbs    number of primitives for each gaussian basis function for each element
// ngbs        number of gaussian basis functions for each element
// scalegbs    scale for each gaussian basis function for each element
// namegbs     name of gaussian basis set for each element
// amgbs       angular momentum for each gaussian basis function for each element
// typgbs      cartesian or spherical gaussian basis type for each element
// coeffgbs    contraction coefficient for each gaussian basis function for each element
// expgbs      primitive exponent for each gaussian basis function for each element

enum class AngMom
{
    S, // = 0
    P, // = 1
    D, // = 2
    F, // = 3
    G, // = 4
    H, // = 5
    I, // = 6
    K, // = 7
};
enum class BasisType
{
    Cartesian,
    Spherical,
};
constexpr int maxnprim = 30;
MDQC_EXTERN MDQCArray<int> nprimgbs;
MDQC_EXTERN MDQCArray<int> ngbs;
MDQC_EXTERN MDQCArray<realq> scalegbs;
MDQC_EXTERN MDQCArray<std::string> namegbs;
MDQC_EXTERN MDQCArray<AngMom> amgbs;
MDQC_EXTERN MDQCArray<BasisType> typgbs;
MDQC_EXTERN MDQCArray2D<realq,maxnprim> coeffgbs;
MDQC_EXTERN MDQCArray2D<realq,maxnprim> expgbs;
}
