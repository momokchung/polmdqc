/////////////////////////////////////////////////////////////
//                                                         //
//  kiprop.h  --  improper dihedral forcefield parameters  //
//                                                         //
/////////////////////////////////////////////////////////////

// maxndi   maximum number of improper dihedral parameter entries
// dcon     force constant parameters for improper dihedrals
// tdi      ideal dihedral angle values for improper dihedrals
// kdi      string of atom classes for improper dihedral angles


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxndi;
QCMD_EXTERN std::vector<double> dcon;
QCMD_EXTERN std::vector<double> tdi;
QCMD_EXTERN std::vector<std::string> kdi;
