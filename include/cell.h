///////////////////////////////////////////////////////
//                                                   //
//  cell.h  --  replicated cell periodic boundaries  //
//                                                   //
///////////////////////////////////////////////////////

// ncell    total number of cell replicates for periodic boundaries
// icell    offset along axes for each replicate periodic cell
// xcell    length of the a-axis of the complete replicated cell
// ycell    length of the b-axis of the complete replicated cell
// zcell    length of the c-axis of the complete replicated cell
// xcell2   half the length of the a-axis of the replicated cell
// ycell2   half the length of the b-axis of the replicated cell
// zcell2   half the length of the c-axis of the replicated cell


#pragma once
#include "macro.h"

QCMD_EXTERN int ncell;
QCMD_EXTERN std::vector<std::vector<int>> icell;
QCMD_EXTERN double xcell;
QCMD_EXTERN double ycell;
QCMD_EXTERN double zcell;
QCMD_EXTERN double xcell2;
QCMD_EXTERN double ycell2;
QCMD_EXTERN double zcell2;
