// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  mutant  --  free energy calculation hybrid atoms  //
//                                                    //
////////////////////////////////////////////////////////

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

MDQC_EXTERN int nmut;
MDQC_EXTERN int vcouple;
MDQC_EXTERN MDQCArray<int> imut;
MDQC_EXTERN MDQCArray<int> type0;
MDQC_EXTERN MDQCArray<int> class0;
MDQC_EXTERN MDQCArray<int> type1;
MDQC_EXTERN MDQCArray<int> class1;
MDQC_EXTERN real lambda;
MDQC_EXTERN real vlambda;
MDQC_EXTERN real elambda;
MDQC_EXTERN real tlambda;
MDQC_EXTERN real scexp;
MDQC_EXTERN real scalpha;
MDQC_EXTERN MDQCArray<bool> mut;
}
