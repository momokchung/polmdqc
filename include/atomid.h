/////////////////////////////////////////////////////////
//                                                     //
//  atomid.h  --  atomic properties for current atoms  //
//                                                     //
/////////////////////////////////////////////////////////

// tag       integer atom labels from input coordinates file
// class     atom class number for each atom in the system
// atomic    atomic number for each atom in the system
// valence   valence number for each atom in the system
// mass      atomic weight for each atom in the system
// name      atom name for each atom in the system
// story     descriptive type for each atom in system


#pragma once
#include "macro.h"
#include "sizes.h"
#include <string>

QCMD_EXTERN int tag[maxatm];
QCMD_EXTERN int atomClass[maxatm];
QCMD_EXTERN int atomic[maxatm];
QCMD_EXTERN int valence[maxatm];
QCMD_EXTERN double mass[maxatm];
QCMD_EXTERN std::string name[maxatm];
QCMD_EXTERN std::string story[maxatm];
