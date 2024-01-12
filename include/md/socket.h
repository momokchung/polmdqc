// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  socket  --  socket communication control parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// skttyp      socket information type (1=DYN, 2=OPT)
// cstep       current dynamics or optimization step number
// cdt         current dynamics cumulative simulation time
// cenergy     current potential energy from simulation
// sktstart    logical flag to indicate socket initialization
// sktstop     logical flag to indicate socket shutdown
// use_socket  logical flag governing use of external sockets

MDQC_EXTERN int skttyp;
MDQC_EXTERN int cstep;
MDQC_EXTERN real cdt;
MDQC_EXTERN real cenergy;
MDQC_EXTERN bool sktstart;
MDQC_EXTERN bool sktstop;
MDQC_EXTERN bool use_socket;
}
