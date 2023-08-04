///////////////////////////////////////////////////////////
//                                                       //
//  mpole.h  --  atomic multipoles in current structure  //
//                                                       //
///////////////////////////////////////////////////////////

// maxpole   max components (monopole=1,dipole=4,quadrupole=13)
// npole     total number of multipole sites in the system
// ipole     number of the atom for each multipole site
// polsiz    number of multipole components for each atom
// pollist   multipole site for each atom (-1=no multipole)
// zaxis     number of the z-axis defining atom for each atom (starting index=1)
// xaxis     number of the x-axis defining atom for each atom (starting index=1)
// yaxis     number of the y-axis defining atom for each atom (starting index=1)
// pole      local frame Cartesian multipoles for each atom
// rpole     global frame Cartesian multipoles for each atom
// mono0     original atomic monopole values for charge flux
// polaxe    local coordinate frame type for each atom


#pragma once
#include "macro.h"
#include <string>
#include <vector>

constexpr int maxpole = 13;
QCMD_EXTERN int npole;
QCMD_EXTERN std::vector<int> ipole;
QCMD_EXTERN std::vector<int> polsiz;
QCMD_EXTERN std::vector<int> pollist;
QCMD_EXTERN std::vector<int> zaxis;
QCMD_EXTERN std::vector<int> xaxis;
QCMD_EXTERN std::vector<int> yaxis;
QCMD_EXTERN std::vector<std::vector<double>> pole;
QCMD_EXTERN std::vector<std::vector<double>> rpole;
QCMD_EXTERN std::vector<double>  mono0;
QCMD_EXTERN std::vector<std::string> polaxe;
