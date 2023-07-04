/////////////////////////////////////////////////////////
//                                                     //
//  params.h  --  force field parameter file contents  //
//                                                     //
/////////////////////////////////////////////////////////

// maxprm    maximum number of lines in the parameter file
// nprm      number of nonblank lines in the parameter file
// prmline   contents of each individual parameter file line


#pragma once
#include "macro.h"
#include <string>

const int maxprm = 25000;
QCMD_EXTERN int nprm;
QCMD_EXTERN std::string prmline[maxprm];
