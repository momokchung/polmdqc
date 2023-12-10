// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  inform  --  program I/O and flow control values  //
//                                                   //
///////////////////////////////////////////////////////

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

constexpr int maxask = 5;
MDQC_EXTERN int gpucard,digits;
MDQC_EXTERN int iprint,iwrite;
MDQC_EXTERN int isend;
MDQC_EXTERN bool verbose,debug;
MDQC_EXTERN bool silent,holdup;
MDQC_EXTERN bool informAbort;
}
