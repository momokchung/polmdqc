// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  kmulti  --  atomic multipole forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// maxnmp   maximum number of atomic multipole parameter entries
// multip   atomic monopole, dipole and quadrupole values
// mpaxis   type of local axis definition for atomic multipoles
// kmp      string of atom types for atomic multipoles

MDQC_EXTERN int maxnmp;
MDQC_EXTERN std::vector<std::vector<double>> multip;
MDQC_EXTERN std::vector<std::string> mpaxis;
MDQC_EXTERN std::vector<std::string> kmp;
}
