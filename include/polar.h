//////////////////////////////////////////////////////////
//                                                      //
//  polar.h  --  polarization & induced dipole moments  //
//                                                      //
//////////////////////////////////////////////////////////

// npolar    total number of polarizable sites in the system
// ipolar    number of the atom for each polarizable site
// jpolar    index into polarization parameter matrix for each atom
// polarity  dipole polarizability for each atom site (Ang**3)
// thole     Thole polarization damping value for each atom
// tholed    Thole direct polarization damping value for each atom
// pdamp     value of polarizability scale factor for each atom
// thlval    Thole damping parameter value for each atom type pair
// thdval    alternate Thole direct damping value for atom type pair
// udir      direct induced dipole components for each atom site
// udirp     direct induced dipoles in field used for energy terms
// udirs     direct GK or PB induced dipoles for each atom site
// udirps    direct induced dipoles in field used for GK or PB energy
// uind      mutual induced dipole components for each atom site
// uinp      mutual induced dipoles in field used for energy terms
// uinds     mutual GK or PB induced dipoles for each atom site
// uinps     mutual induced dipoles in field used for GK or PB energy
// uexact    exact SCF induced dipoles to full numerical precision
// douind    flag to allow induced dipoles at each atom site


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int npolar;
QCMD_EXTERN std::vector<int> ipolar;
QCMD_EXTERN std::vector<int> jpolar;
QCMD_EXTERN std::vector<double> polarity;
QCMD_EXTERN std::vector<double> thole;
QCMD_EXTERN std::vector<double> tholed;
QCMD_EXTERN std::vector<double> pdamp;
QCMD_EXTERN std::vector<std::vector<double>> thlval;
QCMD_EXTERN std::vector<std::vector<double>> thdval;
QCMD_EXTERN std::vector<std::vector<double>> udir;
QCMD_EXTERN std::vector<std::vector<double>> udirp;
QCMD_EXTERN std::vector<std::vector<double>> udirs;
QCMD_EXTERN std::vector<std::vector<double>> udirps;
QCMD_EXTERN std::vector<std::vector<double>> uind;
QCMD_EXTERN std::vector<std::vector<double>> uinp;
QCMD_EXTERN std::vector<std::vector<double>> uinds;
QCMD_EXTERN std::vector<std::vector<double>> uinps;
QCMD_EXTERN std::vector<std::vector<double>> uexact;
QCMD_EXTERN std::vector<bool> douind;
