/////////////////////////////////////////////////////////////
//                                                         //
//  kvdws.h  --  van der Waals term forcefield parameters  //
//                                                         //
/////////////////////////////////////////////////////////////

// rad      van der Waals radius parameter for each atom class
// eps      van der Waals well depth parameter for each atom class
// rad4     van der Waals radius parameter in 1-4 interactions
// eps4     van der Waals well depth parameter in 1-4 interactions
// reduct   van der Waals reduction factor for each atom class


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<double> rad;
QCMD_EXTERN std::vector<double> eps;
QCMD_EXTERN std::vector<double> rad4;
QCMD_EXTERN std::vector<double> eps4;
QCMD_EXTERN std::vector<double> reduct;
