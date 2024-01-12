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

MDQC_EXTERN int tag[maxatm];
MDQC_EXTERN int atomClass[maxatm];
MDQC_EXTERN int atomic[maxatm];
MDQC_EXTERN int valence[maxatm];
MDQC_EXTERN real mass[maxatm];
MDQC_EXTERN std::string name[maxatm];
MDQC_EXTERN std::string story[maxatm];
}
