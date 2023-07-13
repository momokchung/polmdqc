//////////////////////////////////////////////////////////
//                                                      //
//  khbond.h  --  H-bonding term forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// maxnhb   maximum number of hydrogen bonding pair entries
// radhb    radius parameter for hydrogen bonding pairs
// epshb    well depth parameter for hydrogen bonding pairs
// khb      string of atom types for hydrogen bonding pairs


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnhb;
QCMD_EXTERN std::vector<double> radhb;
QCMD_EXTERN std::vector<double> epshb;
QCMD_EXTERN std::vector<std::string> khb;
