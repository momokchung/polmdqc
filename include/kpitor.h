/////////////////////////////////////////////////////////////
//                                                         //
//  kpitor.h  --  pi-system torsion forcefield parameters  //
//                                                         //
/////////////////////////////////////////////////////////////

// maxnpt   maximum number of pi-system torsion parameter entries
// ptcon    force constant parameters for pi-system torsions
// kpt      string of atom classes for pi-system torsion terms


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnpt;
QCMD_EXTERN std::vector<double> ptcon;
QCMD_EXTERN std::vector<std::string> kpt;
