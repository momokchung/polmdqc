// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  virial  --  components of internal virial tensor  //
//                                                    //
////////////////////////////////////////////////////////

// vir         total internal virial Cartesian tensor components
// use_virial  logical flag governing use of virial computation

MDQC_EXTERN real vir[3][3];
MDQC_EXTERN bool use_virial;
}
