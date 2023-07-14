////////////////////////////////////////////////////////////
//                                                        //
//  kvdwpr.h  --  special pair vdw forcefield parameters  //
//                                                        //
////////////////////////////////////////////////////////////

// maxnvp   maximum number of special pair van der Waals entries
// radpr    radius parameter for special van der Waals pairs
// epspr    well depth parameter for special van der Waals pairs
// kvpr     string of atom classes for special van der Waals pairs


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnvp;
QCMD_EXTERN std::vector<double> radpr;
QCMD_EXTERN std::vector<double> epspr;
QCMD_EXTERN std::vector<std::string> kvpr;
