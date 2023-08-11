////////////////////////////////////////////////////////////
//                                                        //
//  files.h  --  name & number of current structure file  //
//                                                        //
////////////////////////////////////////////////////////////

// nprior     number of previously existing cycle files
// ldir       length in characters of the directory name
// leng       length in characters of the base filename
// filename   full filename including any extension or version
// outfile    output filename used for intermediate results


#pragma once
#include "macro.h"
#include <fstream>
#include <string>

QCMD_EXTERN int nprior;
QCMD_EXTERN int ldir,leng;
QCMD_EXTERN std::string filename;
QCMD_EXTERN std::string outfile;
QCMD_EXTERN std::ifstream ffile;
