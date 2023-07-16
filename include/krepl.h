//////////////////////////////////////////////////////////
//                                                      //
//  krepl.h  --  Pauli repulsion forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// prsiz     Pauli repulsion size value for each atom class
// prdmp     alpha Pauli repulsion parameter for each atom class
// prele     number of valence electrons for each atom class


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<double> prsiz;
QCMD_EXTERN std::vector<double> prdmp;
QCMD_EXTERN std::vector<double> prele;
