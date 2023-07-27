///////////////////////////////////////////////////////
//                                                   //
//  kdipol.h  --  bond dipole forcefield parameters  //
//                                                   //
///////////////////////////////////////////////////////

// maxnd    maximum number of bond dipole parameter entries
// maxnd5   maximum number of 5-membered ring dipole entries
// maxnd4   maximum number of 4-membered ring dipole entries
// maxnd3   maximum number of 3-membered ring dipole entries
// dpl      dipole moment parameters for bond dipoles
// dpl5     dipole moment parameters for 5-ring dipoles
// dpl4     dipole moment parameters for 4-ring dipoles
// dpl3     dipole moment parameters for 3-ring dipoles
// pos      dipole position parameters for bond dipoles
// pos5     dipole position parameters for 5-ring dipoles
// pos4     dipole position parameters for 4-ring dipoles
// pos3     dipole position parameters for 3-ring dipoles
// kd       string of atom classes for bond dipoles
// kd5      string of atom classes for 5-ring dipoles
// kd4      string of atom classes for 4-ring dipoles
// kd3      string of atom classes for 3-ring dipoles


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnd;
QCMD_EXTERN int maxnd5;
QCMD_EXTERN int maxnd4;
QCMD_EXTERN int maxnd3;
QCMD_EXTERN std::vector<double> dpl;
QCMD_EXTERN std::vector<double> dpl5;
QCMD_EXTERN std::vector<double> dpl4;
QCMD_EXTERN std::vector<double> dpl3;
QCMD_EXTERN std::vector<double> pos;
QCMD_EXTERN std::vector<double> pos5;
QCMD_EXTERN std::vector<double> pos4;
QCMD_EXTERN std::vector<double> pos3;
QCMD_EXTERN std::vector<std::string> kd;
QCMD_EXTERN std::vector<std::string> kd5;
QCMD_EXTERN std::vector<std::string> kd4;
QCMD_EXTERN std::vector<std::string> kd3;