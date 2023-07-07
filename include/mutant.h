//////////////////////////////////////////////////////////
//                                                      //
//  mutant.h  --  free energy calculation hybrid atoms  //
//                                                      //
//////////////////////////////////////////////////////////

// nmut       number of atoms mutated from initial to final state
// vcouple    van der Waals lambda type (0=decouple, 1=annihilate)
// imut       atom sites differing in initial and final state
// type0      atom type of each atom in the initial state system
// class0     atom class of each atom in the initial state system
// type1      atom type of each atom in the final state system
// class1     atom class of each atom in the final state system
// lambda     generic weighting between initial and final states
// vlambda    state weighting value for van der Waals potentials
// elambda    state weighting value for electrostatic potentials
// tlambda    state weighting value for torsional potential
// scexp      scale factor for soft core buffered 14-7 potential
// scalpha    scale factor for soft core buffered 14-7 potential
// mut        true if an atom is to be mutated, false otherwise


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nmut;
QCMD_EXTERN int vcouple;
QCMD_EXTERN std::vector<int> imut;
QCMD_EXTERN std::vector<int> type0;
QCMD_EXTERN std::vector<int> class0;
QCMD_EXTERN std::vector<int> type1;
QCMD_EXTERN std::vector<int> class1;
QCMD_EXTERN double lambda;
QCMD_EXTERN double vlambda;
QCMD_EXTERN double elambda;
QCMD_EXTERN double tlambda;
QCMD_EXTERN double scexp;
QCMD_EXTERN double scalpha;
QCMD_EXTERN std::vector<bool> mut;
