//////////////////////////////////////////////////////////
//                                                      //
//  chgtrn.h  --  charge transfer in current structure  //
//                                                      //
//////////////////////////////////////////////////////////

// nct       total number of dispersion sites in the system
// chgct     charge for charge transfer at each multipole site
// dmpct     charge transfer damping factor at each multipole site


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nct;
QCMD_EXTERN std::vector<double> chgct;
QCMD_EXTERN std::vector<double> dmpct;
