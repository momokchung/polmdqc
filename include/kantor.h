/////////////////////////////////////////////////////////
//                                                     //
//  kantor.h  --  angle-torsion forcefield parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// maxnat   maximum number of angle-torsion parameter entries
// atcon    torsional amplitude parameters for angle-torsion
// kat      string of atom classes for angle-torsion terms


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnat;
QCMD_EXTERN std::vector<std::vector<double>> atcon;
QCMD_EXTERN std::vector<std::string> kat;
