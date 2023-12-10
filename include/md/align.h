// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  align  --  information for structure superposition  //
//                                                      //
//////////////////////////////////////////////////////////

// nfit    number of atoms to use in superimposing two structures
// ifit    atom numbers of pairs of atoms to be superimposed
// wfit    weights assigned to atom pairs during superposition

MDQC_EXTERN int nfit;
MDQC_EXTERN std::vector<std::vector<int>> ifit;
MDQC_EXTERN std::vector<double> wfit;
}
