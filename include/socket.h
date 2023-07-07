/////////////////////////////////////////////////////////////
//                                                         //
//  socket.h  --  socket communication control parameters  //
//                                                         //
/////////////////////////////////////////////////////////////

// skttyp      socket information type (1=DYN, 2=OPT)
// cstep       current dynamics or optimization step number
// cdt         current dynamics cumulative simulation time
// cenergy     current potential energy from simulation
// sktstart    logical flag to indicate socket initialization
// sktstop     logical flag to indicate socket shutdown
// use_socket  logical flag governing use of external sockets


#pragma once
#include "macro.h"

QCMD_EXTERN int skttyp;
QCMD_EXTERN int cstep;
QCMD_EXTERN double cdt;
QCMD_EXTERN double cenergy;
QCMD_EXTERN bool sktstart;
QCMD_EXTERN bool sktstop;
QCMD_EXTERN bool use_socket;
