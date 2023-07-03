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
#include <string>

const int maxask = 5;
extern int gpucard,digits;
extern int iprint,iwrite;
extern int isend;
extern bool verbose,debug;
extern bool silent,holdup;
extern bool informAbort;
