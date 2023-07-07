/////////////////////////////////////////////////////////
//                                                     //
//  inform.h  --  program I/O and flow control values  //
//                                                     //
/////////////////////////////////////////////////////////

// maxask    maximum number of queries for interactive input
// gpucard   integer flag for GPU use (0=no GPU, 1=GPU present)
// digits    decimal places output for energy and coordinates
// iprint    steps between status printing (0=no printing)
// iwrite    steps between coordinate saves (0=no saves)
// isend     steps between socket communication (0=no sockets)
// verbose   logical flag to turn on extra information printing
// debug     logical flag to turn on extensive debug printing
// silent    logical flag to turn off all information printing
// holdup    logical flag to wait for carriage return on exit
// abort     logical flag to stop execution at next chance


#pragma once
#include "macro.h"
#include <string>

constexpr int maxask = 5;
QCMD_EXTERN int gpucard,digits;
QCMD_EXTERN int iprint,iwrite;
QCMD_EXTERN int isend;
QCMD_EXTERN bool verbose,debug;
QCMD_EXTERN bool silent,holdup;
QCMD_EXTERN bool informAbort;
