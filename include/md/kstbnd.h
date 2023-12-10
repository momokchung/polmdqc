// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  kstbnd  --  stretch-bend forcefield parameters  //
//                                                  //
//////////////////////////////////////////////////////

// maxnsb   maximum number of stretch-bend parameter entries
// stbn     force constant parameters for stretch-bend terms
// ksb      string of atom classes for stretch-bend terms

MDQC_EXTERN int maxnsb;
MDQC_EXTERN std::vector<std::vector<double>> stbn;
MDQC_EXTERN std::vector<std::string> ksb;
}
