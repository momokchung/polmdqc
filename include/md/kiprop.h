// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  kiprop  --  improper dihedral forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// maxndi   maximum number of improper dihedral parameter entries
// dcon     force constant parameters for improper dihedrals
// tdi      ideal dihedral angle values for improper dihedrals
// kdi      string of atom classes for improper dihedral angles

MDQC_EXTERN int maxndi;
MDQC_EXTERN std::vector<real> dcon;
MDQC_EXTERN std::vector<real> tdi;
MDQC_EXTERN std::vector<std::string> kdi;
}
