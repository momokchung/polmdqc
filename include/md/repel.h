// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "mpole.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  repel  --  Pauli repulsion for current structure  //
//                                                    //
////////////////////////////////////////////////////////

// nrep      total number of repulsion sites in the system
// irep      number of the atom for each repulsion site
// replist   repulsion multipole site for each atom (-1=none)
// sizpr     Pauli repulsion size parameter value for each atom
// dmppr     Pauli repulsion alpha damping value for each atom
// elepr     Pauli repulsion valence electrons for each atom
// repole    repulsion Cartesian multipoles in the local frame
// rrepole   repulsion Cartesian multipoles in the global frame

MDQC_EXTERN int nrep;
MDQC_EXTERN MDQCArray<int> irep;
MDQC_EXTERN MDQCArray<int> replist;
MDQC_EXTERN MDQCArray<real> sizpr;
MDQC_EXTERN MDQCArray<real> dmppr;
MDQC_EXTERN MDQCArray<real> elepr;
MDQC_EXTERN MDQCArray2D<real,maxpole> repole;
MDQC_EXTERN MDQCArray2D<real,maxpole> rrepole;
}
