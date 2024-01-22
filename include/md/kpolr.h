// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "sizes.h"

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

MDQC_EXTERN MDQCArray2D<int,maxval> pgrp;
MDQC_EXTERN MDQCArray<real> polr;
MDQC_EXTERN MDQCArray<real> athl;
MDQC_EXTERN MDQCArray<real> dthl;
}
