// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  kdipol  --  bond dipole forcefield parameters  //
//                                                 //
/////////////////////////////////////////////////////

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

MDQC_EXTERN int maxnd;
MDQC_EXTERN int maxnd5;
MDQC_EXTERN int maxnd4;
MDQC_EXTERN int maxnd3;
MDQC_EXTERN MDQCArray<real> dpl;
MDQC_EXTERN MDQCArray<real> dpl5;
MDQC_EXTERN MDQCArray<real> dpl4;
MDQC_EXTERN MDQCArray<real> dpl3;
MDQC_EXTERN MDQCArray<real> pos;
MDQC_EXTERN MDQCArray<real> pos5;
MDQC_EXTERN MDQCArray<real> pos4;
MDQC_EXTERN MDQCArray<real> pos3;
MDQC_EXTERN MDQCArray<std::string> kd;
MDQC_EXTERN MDQCArray<std::string> kd5;
MDQC_EXTERN MDQCArray<std::string> kd4;
MDQC_EXTERN MDQCArray<std::string> kd3;
}
