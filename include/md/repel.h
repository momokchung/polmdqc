// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  repel  --  Pauli repulsion for current structure  //
//                                                    //
////////////////////////////////////////////////////////

// nrep      total number of repulsion sites in the system
// irep      number of the atom for each repulsion site
// replist   repulsion multipole site for each atom (0=none)
// sizpr     Pauli repulsion size parameter value for each atom
// dmppr     Pauli repulsion alpha damping value for each atom
// elepr     Pauli repulsion valence electrons for each atom
// repole    repulsion Cartesian multipoles in the local frame
// rrepole   repulsion Cartesian multipoles in the global frame

MDQC_EXTERN int nrep;
MDQC_EXTERN std::vector<int> irep;
MDQC_EXTERN std::vector<int> replist;
MDQC_EXTERN std::vector<double> sizpr;
MDQC_EXTERN std::vector<double> dmppr;
MDQC_EXTERN std::vector<double> elepr;
MDQC_EXTERN std::vector<std::vector<double>> repole;
MDQC_EXTERN std::vector<std::vector<double>> rrepole;
}
