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
#include <string>

const int maxarg = 20;
extern int narg;
extern bool listarg[maxarg];
extern std::string arg[maxarg];
