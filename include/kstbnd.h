////////////////////////////////////////////////////////
//                                                    //
//  kstbnd.h  --  stretch-bend forcefield parameters  //
//                                                    //
////////////////////////////////////////////////////////

// maxnsb   maximum number of stretch-bend parameter entries
// stbn     force constant parameters for stretch-bend terms
// ksb      string of atom classes for stretch-bend terms


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnsb;
QCMD_EXTERN std::vector<std::vector<double>> stbn;
QCMD_EXTERN std::vector<std::string> ksb;
