// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "sizes.h"

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  couple  --  atom neighbor connectivity lists  //
//                                                //
////////////////////////////////////////////////////

// n12      number of atoms directly bonded to each atom
// n13      number of atoms in a 1-3 relation to each atom
// n14      number of atoms in a 1-4 relation to each atom
// n15      number of atoms in a 1-5 relation to each atom
// i12      atom numbers of atoms 1-2 connected to each atom
// i13      atom numbers of atoms 1-3 connected to each atom
// i14      atom numbers of atoms 1-4 connected to each atom
// i15      atom numbers of atoms 1-5 connected to each atom

MDQC_EXTERN MDQCArray<int> n12;
MDQC_EXTERN MDQCArray<int> n13;
MDQC_EXTERN MDQCArray<int> n14;
MDQC_EXTERN MDQCArray<int> n15;
MDQC_EXTERN MDQCArray2D<int,maxval> i12;
MDQC_EXTERN MDQCArray2D<int,3*maxval> i13;
MDQC_EXTERN MDQCArray2D<int,9*maxval> i14;
MDQC_EXTERN MDQCArray2D<int,27*maxval> i15;
}
