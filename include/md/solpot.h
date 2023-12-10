// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  solpot  --  solvation term functional form details  //
//                                                      //
//////////////////////////////////////////////////////////

// solvtyp   type of continuum solvation energy model in use
// borntyp   method to be used for the Born radius computation

MDQC_EXTERN std::string solvtyp;
MDQC_EXTERN std::string borntyp;
}
