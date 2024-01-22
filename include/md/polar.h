// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

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
MDQC_EXTERN MDQCArray<int> ipolar;
MDQC_EXTERN MDQCArray<int> jpolar;
MDQC_EXTERN MDQCArray<real> polarity;
MDQC_EXTERN MDQCArray<real> thole;
MDQC_EXTERN MDQCArray<real> tholed;
MDQC_EXTERN MDQCArray<real> pdamp;
MDQC_EXTERN MDQCArray<real> thlval;
MDQC_EXTERN MDQCArray<real> thdval;
MDQC_EXTERN MDQCArray2D<real,3> udir;
MDQC_EXTERN MDQCArray2D<real,3> udirp;
MDQC_EXTERN MDQCArray2D<real,3> udirs;
MDQC_EXTERN MDQCArray2D<real,3> udirps;
MDQC_EXTERN MDQCArray2D<real,3> uind;
MDQC_EXTERN MDQCArray2D<real,3> uinp;
MDQC_EXTERN MDQCArray2D<real,3> uinds;
MDQC_EXTERN MDQCArray2D<real,3> uinps;
MDQC_EXTERN MDQCArray2D<real,3> uexact;
MDQC_EXTERN MDQCArray<bool> douind;
}
