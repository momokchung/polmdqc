// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  params  --  force field parameter file contents  //
//                                                   //
///////////////////////////////////////////////////////

// maxprm    maximum number of lines in the parameter file
// nprm      number of nonblank lines in the parameter file
// prmline   contents of each individual parameter file line

constexpr int maxprm = 25000;
MDQC_EXTERN int nprm;
MDQC_EXTERN std::string prmline[maxprm];
}
