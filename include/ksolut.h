//////////////////////////////////////////////////////////
//                                                      //
//  ksolut.h  --  solvation term forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// pbr      Poisson-Boltzmann radius value for each atom type
// csr      ddCOSMO solvation radius value for each atom type
// gkr      Generalized Kirkwood radius value for each atom type


#pragma once
#include "macro.h"
// #include <string>
// #include <vector>

QCMD_EXTERN std::vector<double> pbr;
QCMD_EXTERN std::vector<double> csr;
QCMD_EXTERN std::vector<double> gkr;
