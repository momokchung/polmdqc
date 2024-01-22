// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "sizes.h"
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  atomid  --  atomic properties for current atoms  //
//                                                   //
///////////////////////////////////////////////////////

// tag       integer atom labels from input coordinates file
// class     atom class number for each atom in the system
// atomic    atomic number for each atom in the system
// valence   valence number for each atom in the system
// mass      atomic weight for each atom in the system
// name      atom name for each atom in the system
// story     descriptive type for each atom in system

MDQC_EXTERN MDQCArray<int> tag;
MDQC_EXTERN MDQCArray<int> atomClass;
MDQC_EXTERN MDQCArray<int> atomic;
MDQC_EXTERN MDQCArray<int> valence;
MDQC_EXTERN MDQCArray<real> mass;
MDQC_EXTERN MDQCArray<std::string> name;
MDQC_EXTERN MDQCArray<std::string> story;
}
