/////////////////////////////////////////////////////////
//                                                     //
//  params.h  --  force field parameter file contents  //
//                                                     //
/////////////////////////////////////////////////////////

// maxprm    maximum number of lines in the parameter file
// nprm      number of nonblank lines in the parameter file
// prmline   contents of each individual parameter file line


#pragma once
#include <string>

const int maxprm = 25000;
extern int nprm;
std::string prmline[maxprm];
