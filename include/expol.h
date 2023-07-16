///////////////////////////////////////////////////////////
//                                                       //
//  expol.h  --  exch-polarization in current structure  //
//                                                       //
///////////////////////////////////////////////////////////

// nexpol     total number of exch polarization sites in system
// kpep       exchange polarization spring constant at each site
// prepep     exchange polarization prefactor at each site
// dmppep     exchange polarization damping alpha at each site
// polscale   scale matrix for use in exchange polarization
// polinv     scale matrix inverse for exchange polarization
// lpep       flag to use exchange polarization at each site


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nexpol;
QCMD_EXTERN std::vector<double> kpep;
QCMD_EXTERN std::vector<double> prepep;
QCMD_EXTERN std::vector<double> dmppep;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> polscale;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> polinv;
QCMD_EXTERN std::vector<bool> lpep;
