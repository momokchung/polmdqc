// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  krepl  --  Pauli repulsion forcefield parameters  //
//                                                    //
////////////////////////////////////////////////////////

// prsiz     Pauli repulsion size value for each atom class
// prdmp     alpha Pauli repulsion parameter for each atom class
// prele     number of valence electrons for each atom class

MDQC_EXTERN MDQCArray<real> prsiz;
MDQC_EXTERN MDQCArray<real> prdmp;
MDQC_EXTERN MDQCArray<real> prele;
}
