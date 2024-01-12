// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  kpolr  --  polarizability forcefield parameters  //
//                                                   //
///////////////////////////////////////////////////////

// pgrp   connected types in polarization group of each atom type
// polr   dipole polarizability parameters for each atom type
// athl   Thole polarization damping value for each atom type
// dthl   alternate Thole direct polarization damping values

MDQC_EXTERN std::vector<std::vector<int>> pgrp;
MDQC_EXTERN std::vector<real> polr;
MDQC_EXTERN std::vector<real> athl;
MDQC_EXTERN std::vector<real> dthl;
}
