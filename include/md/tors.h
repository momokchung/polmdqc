// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  tors  --  torsional angles in current structure  //
//                                                   //
///////////////////////////////////////////////////////

// ntors   total number of torsional angles in the system
// itors   numbers of the atoms in each torsional angle
// tors1   1-fold amplitude and phase for each torsional angle
// tors2   2-fold amplitude and phase for each torsional angle
// tors3   3-fold amplitude and phase for each torsional angle
// tors4   4-fold amplitude and phase for each torsional angle
// tors5   5-fold amplitude and phase for each torsional angle
// tors6   6-fold amplitude and phase for each torsional angle

MDQC_EXTERN int ntors;
MDQC_EXTERN std::vector<std::vector<int>>  itors;
MDQC_EXTERN std::vector<std::vector<double>> tors1;
MDQC_EXTERN std::vector<std::vector<double>> tors2;
MDQC_EXTERN std::vector<std::vector<double>> tors3;
MDQC_EXTERN std::vector<std::vector<double>> tors4;
MDQC_EXTERN std::vector<std::vector<double>> tors5;
MDQC_EXTERN std::vector<std::vector<double>> tors6;
}
