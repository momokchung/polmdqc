////////////////////////////////////////////////////////
//                                                    //
//  keys.h  --  contents of the keyword control file  //
//                                                    //
////////////////////////////////////////////////////////

// maxkey    maximum number of lines in the keyword file
// nkey      number of nonblank lines in the keyword file
// keyline   contents of each individual keyword file line


#pragma once
#include <string>

const int maxkey = 25000;
extern int nkey;
std::string keyline[maxkey];
