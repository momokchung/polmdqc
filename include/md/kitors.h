// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  kitors  --  improper torsion forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// maxnti   maximum number of improper torsion parameter entries
// ti1      torsional parameters for improper 1-fold rotation
// ti2      torsional parameters for improper 2-fold rotation
// ti3      torsional parameters for improper 3-fold rotation
// kti      string of atom classes for improper torsional parameters

MDQC_EXTERN int maxnti;
MDQC_EXTERN std::vector<std::vector<real>> ti1;
MDQC_EXTERN std::vector<std::vector<real>> ti2;
MDQC_EXTERN std::vector<std::vector<real>> ti3;
MDQC_EXTERN std::vector<std::string> kti;
}
