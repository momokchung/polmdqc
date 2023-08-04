//////////////////////////////////////////////////////////
//                                                      //
//  repel.h  --  Pauli repulsion for current structure  //
//                                                      //
//////////////////////////////////////////////////////////

// nrep      total number of repulsion sites in the system
// irep      number of the atom for each repulsion site
// replist   repulsion multipole site for each atom (0=none)
// sizpr     Pauli repulsion size parameter value for each atom
// dmppr     Pauli repulsion alpha damping value for each atom
// elepr     Pauli repulsion valence electrons for each atom
// repole    repulsion Cartesian multipoles in the local frame
// rrepole   repulsion Cartesian multipoles in the global frame


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nrep;
QCMD_EXTERN std::vector<int> irep;
QCMD_EXTERN std::vector<int> replist;
QCMD_EXTERN std::vector<double> sizpr;
QCMD_EXTERN std::vector<double> dmppr;
QCMD_EXTERN std::vector<double> elepr;
QCMD_EXTERN std::vector<std::vector<double>> repole;
QCMD_EXTERN std::vector<std::vector<double>> rrepole;
