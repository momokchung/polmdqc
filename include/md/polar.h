// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  polar  --  polarization & induced dipole moments  //
//                                                    //
////////////////////////////////////////////////////////

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

MDQC_EXTERN int npolar;
MDQC_EXTERN std::vector<int> ipolar;
MDQC_EXTERN std::vector<int> jpolar;
MDQC_EXTERN std::vector<real> polarity;
MDQC_EXTERN std::vector<real> thole;
MDQC_EXTERN std::vector<real> tholed;
MDQC_EXTERN std::vector<real> pdamp;
MDQC_EXTERN std::vector<std::vector<real>> thlval;
MDQC_EXTERN std::vector<std::vector<real>> thdval;
MDQC_EXTERN std::vector<std::vector<real>> udir;
MDQC_EXTERN std::vector<std::vector<real>> udirp;
MDQC_EXTERN std::vector<std::vector<real>> udirs;
MDQC_EXTERN std::vector<std::vector<real>> udirps;
MDQC_EXTERN std::vector<std::vector<real>> uind;
MDQC_EXTERN std::vector<std::vector<real>> uinp;
MDQC_EXTERN std::vector<std::vector<real>> uinds;
MDQC_EXTERN std::vector<std::vector<real>> uinps;
MDQC_EXTERN std::vector<std::vector<real>> uexact;
MDQC_EXTERN std::vector<bool> douind;
}
