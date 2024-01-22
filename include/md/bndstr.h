// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  bndstr  --  bond stretches in the current structure  //
//                                                       //
///////////////////////////////////////////////////////////

// nbond   total number of bond stretches in the system
// ibnd    numbers of the atoms in each bond stretch
// bk      bond stretch force constants (kcal/mole/Ang**2)
// bl      ideal bond length values in Angstroms

MDQC_EXTERN int nbond;
MDQC_EXTERN MDQCArray2D<int,2> ibnd;
MDQC_EXTERN MDQCArray<real> bk;
MDQC_EXTERN MDQCArray<real> bl;
}
