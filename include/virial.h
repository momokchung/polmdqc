//////////////////////////////////////////////////////////
//                                                      //
//  virial.h  --  components of internal virial tensor  //
//                                                      //
//////////////////////////////////////////////////////////

// vir         total internal virial Cartesian tensor components
// use_virial  logical flag governing use of virial computation


#pragma once
#include "macro.h"

QCMD_EXTERN double vir[3][3];
QCMD_EXTERN bool use_virial;
