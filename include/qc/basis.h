// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
//////////////////////////////////////////
//                                      //
//  basis  --  basis set file contents  //
//                                      //
//////////////////////////////////////////

// maxbss    maximum number of lines in the basis file
// nbss      number of nonblank lines in the basis file
// bssline   contents of each individual basis file line

constexpr int maxbss = 25000;
MDQC_EXTERN int nbss;
MDQC_EXTERN std::string bssname;
MDQC_EXTERN std::string bssline[maxbss];
}
