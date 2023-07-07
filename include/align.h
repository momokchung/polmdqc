////////////////////////////////////////////////////////////
//                                                        //
//  align.h  --  information for structure superposition  //
//                                                        //
////////////////////////////////////////////////////////////

// nfit    number of atoms to use in superimposing two structures
// ifit    atom numbers of pairs of atoms to be superimposed
// wfit    weights assigned to atom pairs during superposition


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nfit;
QCMD_EXTERN std::vector<std::vector<int>> ifit;
QCMD_EXTERN std::vector<double> wfit;
