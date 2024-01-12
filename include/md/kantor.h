// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  kantor  --  angle-torsion forcefield parameters  //
//                                                   //
///////////////////////////////////////////////////////

// maxnat   maximum number of angle-torsion parameter entries
// atcon    torsional amplitude parameters for angle-torsion
// kat      string of atom classes for angle-torsion terms

MDQC_EXTERN int maxnat;
MDQC_EXTERN std::vector<std::vector<real>> atcon;
MDQC_EXTERN std::vector<std::string> kat;
}
