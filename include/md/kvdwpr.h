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
//  kvdwpr  --  special pair vdw forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// maxnvp   maximum number of special pair van der Waals entries
// radpr    radius parameter for special van der Waals pairs
// epspr    well depth parameter for special van der Waals pairs
// kvpr     string of atom classes for special van der Waals pairs

MDQC_EXTERN int maxnvp;
MDQC_EXTERN std::vector<real> radpr;
MDQC_EXTERN std::vector<real> epspr;
MDQC_EXTERN std::vector<std::string> kvpr;
}
