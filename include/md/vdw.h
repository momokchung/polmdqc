// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  vdw  --  van der Waals terms in current structure  //
//                                                     //
/////////////////////////////////////////////////////////

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

MDQC_EXTERN int nvdw;
MDQC_EXTERN std::vector<int> ivdw;
MDQC_EXTERN std::vector<int> jvdw;
MDQC_EXTERN std::vector<int> mvdw;
MDQC_EXTERN std::vector<int> ired;
MDQC_EXTERN std::vector<real> kred;
MDQC_EXTERN std::vector<real> xred;
MDQC_EXTERN std::vector<real> yred;
MDQC_EXTERN std::vector<real> zred;
MDQC_EXTERN std::vector<std::vector<real>> radmin;
MDQC_EXTERN std::vector<std::vector<real>> epsilon;
MDQC_EXTERN std::vector<std::vector<real>> radmin4;
MDQC_EXTERN std::vector<std::vector<real>> epsilon4;
MDQC_EXTERN std::vector<std::vector<real>> radhbnd;
MDQC_EXTERN std::vector<std::vector<real>> epshbnd;
}
