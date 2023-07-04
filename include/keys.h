////////////////////////////////////////////////////////
//                                                    //
//  keys.h  --  contents of the keyword control file  //
//                                                    //
////////////////////////////////////////////////////////

// maxkey    maximum number of lines in the keyword file
// nkey      number of nonblank lines in the keyword file
// keyline   contents of each individual keyword file line


#pragma once
#include "macro.h"
#include <string>

const int maxkey = 25000;
QCMD_EXTERN int nkey;
QCMD_EXTERN std::string keyline[maxkey];
