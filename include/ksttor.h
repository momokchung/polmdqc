///////////////////////////////////////////////////////////
//                                                       //
//  ksttor.h  --  stretch-torsion forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// maxnbt   maximum number of stretch-torsion parameter entries
// btcon    torsional amplitude parameters for stretch-torsion
// kbt      string of atom classes for stretch-torsion terms


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnbt;
QCMD_EXTERN std::vector<std::vector<double>> btcon;
QCMD_EXTERN std::vector<std::string> kbt;
