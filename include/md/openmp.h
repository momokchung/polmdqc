// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  openmp  --  OpenMP processor and thread values  //
//                                                  //
//////////////////////////////////////////////////////

// nproc     number of processors available to OpenMP
// nthread   number of threads to be used with OpenMP 

MDQC_EXTERN int nproc;
MDQC_EXTERN int nthread;
}
