///////////////////////////////////////////////////////////
//                                                       //
//  vdw.h  --  van der Waals terms in current structure  //
//                                                       //
///////////////////////////////////////////////////////////

// nvdw       total number van der Waals sites in the system
// ivdw       number of the atom for each van der Waals site
// jvdw       index into the vdw parameter matrix for each atom
// mvdw       index into the vdw parameter matrix for each class
// ired       attached atom from which reduction factor is applied
// kred       value of reduction factor parameter for each atom
// xred       reduced x-coordinate for each atom in the system
// yred       reduced y-coordinate for each atom in the system
// zred       reduced z-coordinate for each atom in the system
// radmin     minimum energy distance for each atom class pair
// epsilon    well depth parameter for each atom class pair
// radmin4    minimum energy distance for 1-4 interaction pairs
// epsilon4   well depth parameter for 1-4 interaction pairs
// radhbnd    minimum energy distance for hydrogen bonding pairs
// epshbnd    well depth parameter for hydrogen bonding pairs


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nvdw;
QCMD_EXTERN std::vector<int> ivdw;
QCMD_EXTERN std::vector<int> jvdw;
QCMD_EXTERN std::vector<int> mvdw;
QCMD_EXTERN std::vector<int> ired;
QCMD_EXTERN std::vector<double> kred;
QCMD_EXTERN std::vector<double> xred;
QCMD_EXTERN std::vector<double> yred;
QCMD_EXTERN std::vector<double> zred;
QCMD_EXTERN std::vector<std::vector<double>> radmin;
QCMD_EXTERN std::vector<std::vector<double>> epsilon;
QCMD_EXTERN std::vector<std::vector<double>> radmin4;
QCMD_EXTERN std::vector<std::vector<double>> epsilon4;
QCMD_EXTERN std::vector<std::vector<double>> radhbnd;
QCMD_EXTERN std::vector<std::vector<double>> epshbnd;
