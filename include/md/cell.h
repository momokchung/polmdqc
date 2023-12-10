// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  cell  --  replicated cell periodic boundaries  //
//                                                 //
/////////////////////////////////////////////////////

// ncell    total number of cell replicates for periodic boundaries
// icell    offset along axes for each replicate periodic cell
// xcell    length of the a-axis of the complete replicated cell
// ycell    length of the b-axis of the complete replicated cell
// zcell    length of the c-axis of the complete replicated cell
// xcell2   half the length of the a-axis of the replicated cell
// ycell2   half the length of the b-axis of the replicated cell
// zcell2   half the length of the c-axis of the replicated cell

MDQC_EXTERN int ncell;
MDQC_EXTERN std::vector<std::vector<int>> icell;
MDQC_EXTERN double xcell;
MDQC_EXTERN double ycell;
MDQC_EXTERN double zcell;
MDQC_EXTERN double xcell2;
MDQC_EXTERN double ycell2;
MDQC_EXTERN double zcell2;
}
