// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <fstream>
#include <string>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  files  --  name & number of current structure file  //
//                                                      //
//////////////////////////////////////////////////////////

// nprior     number of previously existing cycle files
// ldir       length in characters of the directory name
// leng       length in characters of the base filename
// filename   full filename including any extension or version
// outfile    output filename used for intermediate results

MDQC_EXTERN int nprior;
MDQC_EXTERN int ldir,leng;
MDQC_EXTERN std::string filename;
MDQC_EXTERN std::string outfile;
MDQC_EXTERN std::ifstream ffile;
}
