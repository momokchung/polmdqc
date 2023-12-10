// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  ksolut  --  solvation term forcefield parameters  //
//                                                    //
////////////////////////////////////////////////////////

// pbr      Poisson-Boltzmann radius value for each atom type
// csr      ddCOSMO solvation radius value for each atom type
// gkr      Generalized Kirkwood radius value for each atom type

MDQC_EXTERN std::vector<double> pbr;
MDQC_EXTERN std::vector<double> csr;
MDQC_EXTERN std::vector<double> gkr;
}
