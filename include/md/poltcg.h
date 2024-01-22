// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  poltcg  --  induced dipoles via the TCG solver  //
//                                                  //
//////////////////////////////////////////////////////

// tcgorder   total number of TCG iterations to be used
// tcgnab     number of mutual induced dipole components
// tcgpeek    value of acceleration factor for TCG peek step
// uad        left-hand side mutual induced d-dipoles
// uap        left-hand side mutual induced p-dipoles
// ubd        right-hand side mutual induced d-dipoles
// ubp        right-hand side mutual induced p-dipoles
// tcgguess   flag to use initial TCG based on direct field

MDQC_EXTERN int tcgorder;
MDQC_EXTERN int tcgnab;
MDQC_EXTERN real tcgpeek;
MDQC_EXTERN MDQCArray2D<real,3> uad;
MDQC_EXTERN MDQCArray2D<real,3> uap;
MDQC_EXTERN MDQCArray2D<real,3> ubd;
MDQC_EXTERN MDQCArray2D<real,3> ubp;
MDQC_EXTERN bool tcgguess;
}
