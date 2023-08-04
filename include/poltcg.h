////////////////////////////////////////////////////////
//                                                    //
//  poltcg.h  --  induced dipoles via the TCG solver  //
//                                                    //
////////////////////////////////////////////////////////

// tcgorder   total number of TCG iterations to be used
// tcgnab     number of mutual induced dipole components
// tcgpeek    value of acceleration factor for TCG peek step
// uad        left-hand side mutual induced d-dipoles
// uap        left-hand side mutual induced p-dipoles
// ubd        right-hand side mutual induced d-dipoles
// ubp        right-hand side mutual induced p-dipoles
// tcgguess   flag to use initial TCG based on direct field


#pragma once
#include "macro.h"

QCMD_EXTERN int tcgorder;
QCMD_EXTERN int tcgnab;
QCMD_EXTERN double tcgpeek;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> uad;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> uap;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> ubd;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> ubp;
QCMD_EXTERN bool tcgguess;
