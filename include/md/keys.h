// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  keys  --  contents of the keyword control file  //
//                                                  //
//////////////////////////////////////////////////////

// maxkey    maximum number of lines in the keyword file
// nkey      number of nonblank lines in the keyword file
// keyline   contents of each individual keyword file line

constexpr int maxkey = 25000;
MDQC_EXTERN int nkey;
MDQC_EXTERN std::string keyline[maxkey];
}
