// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

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

MDQC_EXTERN MDQCArray<real> pbr;
MDQC_EXTERN MDQCArray<real> csr;
MDQC_EXTERN MDQCArray<real> gkr;
}
