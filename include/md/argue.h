// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  argue  --  command line arguments at run time  //
//                                                 //
/////////////////////////////////////////////////////

//  maxarg    maximum number of command line arguments
//  narg      number of command line arguments to the program
//  listarg   flag to mark available command line arguments
//  arg       strings containing the command line arguments

constexpr int maxarg = 20;
MDQC_EXTERN int narg;
MDQC_EXTERN bool listarg[maxarg];
MDQC_EXTERN std::string arg[maxarg];
}
