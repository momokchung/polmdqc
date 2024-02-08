// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  molcul  --  individual molecules in current system  //
//                                                      //
//////////////////////////////////////////////////////////

// nmol      total number of separate molecules in the system
// imol      first and last atom of each molecule in the list
// kmol      contiguous list of the atoms in each molecule
// molcule   number of the molecule to which each atom belongs
// totmass   total weight of all the molecules in the system
// molmass   molecular weight for each molecule in the system
// molchg    molecular charge for each molecule in the system
// molmult   molecular multiplicity for each molecule in the system

MDQC_EXTERN int nmol;
MDQC_EXTERN MDQCArray2D<int,2> imol;
MDQC_EXTERN MDQCArray<int> kmol;
MDQC_EXTERN MDQCArray<int> molcule;
MDQC_EXTERN MDQCArray<int> molchg;
MDQC_EXTERN MDQCArray<int> molmult;
MDQC_EXTERN real totmass;
MDQC_EXTERN MDQCArray<real> molmass;
}
