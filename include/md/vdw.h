// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

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
MDQC_EXTERN MDQCArray<int> ivdw;
MDQC_EXTERN MDQCArray<int> jvdw;
MDQC_EXTERN MDQCArray<int> mvdw;
MDQC_EXTERN MDQCArray<int> ired;
MDQC_EXTERN MDQCArray<real> kred;
MDQC_EXTERN MDQCArray<real> xred;
MDQC_EXTERN MDQCArray<real> yred;
MDQC_EXTERN MDQCArray<real> zred;
MDQC_EXTERN MDQCArray<real> radmin;
MDQC_EXTERN MDQCArray<real> epsilon;
MDQC_EXTERN MDQCArray<real> radmin4;
MDQC_EXTERN MDQCArray<real> epsilon4;
MDQC_EXTERN MDQCArray<real> radhbnd;
MDQC_EXTERN MDQCArray<real> epshbnd;
}
