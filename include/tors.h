/////////////////////////////////////////////////////////
//                                                     //
//  tors.h  --  torsional angles in current structure  //
//                                                     //
/////////////////////////////////////////////////////////

// ntors   total number of torsional angles in the system
// itors   numbers of the atoms in each torsional angle
// tors1   1-fold amplitude and phase for each torsional angle
// tors2   2-fold amplitude and phase for each torsional angle
// tors3   3-fold amplitude and phase for each torsional angle
// tors4   4-fold amplitude and phase for each torsional angle
// tors5   5-fold amplitude and phase for each torsional angle
// tors6   6-fold amplitude and phase for each torsional angle


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int ntors;
QCMD_EXTERN std::vector<std::vector<int>>  itors;
QCMD_EXTERN std::vector<std::vector<double>> tors1;
QCMD_EXTERN std::vector<std::vector<double>> tors2;
QCMD_EXTERN std::vector<std::vector<double>> tors3;
QCMD_EXTERN std::vector<std::vector<double>> tors4;
QCMD_EXTERN std::vector<std::vector<double>> tors5;
QCMD_EXTERN std::vector<std::vector<double>> tors6;
