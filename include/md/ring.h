// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  ring  --  number and location of ring structures  //
//                                                    //
////////////////////////////////////////////////////////

// nring3   total number of 3-membered rings in the system
// nring4   total number of 4-membered rings in the system
// nring5   total number of 5-membered rings in the system
// nring6   total number of 6-membered rings in the system
// nring7   total number of 7-membered rings in the system
// iring3   numbers of the atoms involved in each 3-ring
// iring4   numbers of the atoms involved in each 4-ring
// iring5   numbers of the atoms involved in each 5-ring
// iring6   numbers of the atoms involved in each 6-ring
// iring7   numbers of the atoms involved in each 7-ring

MDQC_EXTERN int nring3;
MDQC_EXTERN int nring4;
MDQC_EXTERN int nring5;
MDQC_EXTERN int nring6;
MDQC_EXTERN int nring7;
MDQC_EXTERN MDQCArray2D<int,3> iring3;
MDQC_EXTERN MDQCArray2D<int,4> iring4;
MDQC_EXTERN MDQCArray2D<int,5> iring5;
MDQC_EXTERN MDQCArray2D<int,6> iring6;
MDQC_EXTERN MDQCArray2D<int,7> iring7;
}
