///////////////////////////////////////////////////////
//                                                   //
//  argue.h  --  command line arguments at run time  //
//                                                   //
///////////////////////////////////////////////////////

//  maxarg    maximum number of command line arguments
//  narg      number of command line arguments to the program
//  listarg   flag to mark available command line arguments
//  arg       strings containing the command line arguments


#pragma once
#include "macro.h"
#include <string>

constexpr int maxarg = 20;
QCMD_EXTERN int narg;
QCMD_EXTERN bool listarg[maxarg];
QCMD_EXTERN std::string arg[maxarg];
