////////////////////////////////////////////////////////////
//                                                        //
//  kmulti.h  --  atomic multipole forcefield parameters  //
//                                                        //
////////////////////////////////////////////////////////////

// maxnmp   maximum number of atomic multipole parameter entries
// multip   atomic monopole, dipole and quadrupole values
// mpaxis   type of local axis definition for atomic multipoles
// kmp      string of atom types for atomic multipoles


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnmp;
QCMD_EXTERN std::vector<std::vector<double>> multip;
QCMD_EXTERN std::vector<std::string> mpaxis;
QCMD_EXTERN std::vector<std::string> kmp;
