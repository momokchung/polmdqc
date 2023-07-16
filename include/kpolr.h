/////////////////////////////////////////////////////////
//                                                     //
//  kpolr.h  --  polarizability forcefield parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// pgrp   connected types in polarization group of each atom type
// polr   dipole polarizability parameters for each atom type
// athl   Thole polarization damping value for each atom type
// dthl   alternate Thole direct polarization damping values


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<std::vector<int>> pgrp;
QCMD_EXTERN std::vector<double> polr;
QCMD_EXTERN std::vector<double> athl;
QCMD_EXTERN std::vector<double> dthl;
