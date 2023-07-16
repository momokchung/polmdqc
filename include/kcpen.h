/////////////////////////////////////////////////////////////
//                                                         //
//  kcpen.h  --  charge penetration forcefield parameters  //
//                                                         //
/////////////////////////////////////////////////////////////

// cpele     valence electron magnitude for each atom class
// cpalp     alpha charge penetration parameter for each atom class


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<double> cpele;
QCMD_EXTERN std::vector<double> cpalp;
