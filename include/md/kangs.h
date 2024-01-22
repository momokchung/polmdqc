// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  kangs  --  bond angle bend forcefield parameters  //
//                                                    //
////////////////////////////////////////////////////////

// maxna    maximum number of harmonic angle bend parameter entries
// maxna5   maximum number of 5-membered ring angle bend entries
// maxna4   maximum number of 4-membered ring angle bend entries
// maxna3   maximum number of 3-membered ring angle bend entries
// maxnap   maximum number of in-plane angle bend parameter entries
// maxnaf   maximum number of Fourier angle bend parameter entries
// acon     force constant parameters for harmonic angle bends
// acon5    force constant parameters for 5-ring angle bends
// acon4    force constant parameters for 4-ring angle bends
// acon3    force constant parameters for 3-ring angle bends
// aconp    force constant parameters for in-plane angle bends
// aconf    force constant parameters for Fourier angle bends
// ang      bond angle parameters for harmonic angle bends
// ang5     bond angle parameters for 5-ring angle bends
// ang4     bond angle parameters for 4-ring angle bends
// ang3     bond angle parameters for 3-ring angle bends
// angp     bond angle parameters for in-plane angle bends
// angf     phase shift angle and periodicity for Fourier bends
// ka       string of atom classes for harmonic angle bends
// ka5      string of atom classes for 5-ring angle bends
// ka4      string of atom classes for 4-ring angle bends
// ka3      string of atom classes for 3-ring angle bends
// kap      string of atom classes for in-plane angle bends
// kaf      string of atom classes for Fourier angle bends

MDQC_EXTERN int maxna;
MDQC_EXTERN int maxna5;
MDQC_EXTERN int maxna4;
MDQC_EXTERN int maxna3;
MDQC_EXTERN int maxnap;
MDQC_EXTERN int maxnaf;
MDQC_EXTERN MDQCArray<real> acon;
MDQC_EXTERN MDQCArray<real> acon5;
MDQC_EXTERN MDQCArray<real> acon4;
MDQC_EXTERN MDQCArray<real> acon3;
MDQC_EXTERN MDQCArray<real> aconp;
MDQC_EXTERN MDQCArray<real> aconf;
MDQC_EXTERN MDQCArray2D<real,3> ang;
MDQC_EXTERN MDQCArray2D<real,3> ang5;
MDQC_EXTERN MDQCArray2D<real,3> ang4;
MDQC_EXTERN MDQCArray2D<real,3> ang3;
MDQC_EXTERN MDQCArray2D<real,2> angp;
MDQC_EXTERN MDQCArray2D<real,2> angf;
MDQC_EXTERN MDQCArray<std::string> ka;
MDQC_EXTERN MDQCArray<std::string> ka5;
MDQC_EXTERN MDQCArray<std::string> ka4;
MDQC_EXTERN MDQCArray<std::string> ka3;
MDQC_EXTERN MDQCArray<std::string> kap;
MDQC_EXTERN MDQCArray<std::string> kaf;
}
