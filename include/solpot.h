////////////////////////////////////////////////////////////
//                                                        //
//  solpot.h  --  solvation term functional form details  //
//                                                        //
////////////////////////////////////////////////////////////

// solvtyp   type of continuum solvation energy model in use
// borntyp   method to be used for the Born radius computation


#pragma once
#include "macro.h"
#include <string>

QCMD_EXTERN std::string solvtyp;
QCMD_EXTERN std::string borntyp;
