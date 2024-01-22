// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  kbonds  --  bond stretching forcefield parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// maxnb   maximum number of bond stretch parameter entries
// maxnb5  maximum number of 5-membered ring bond stretch entries
// maxnb4  maximum number of 4-membered ring bond stretch entries
// maxnb3  maximum number of 3-membered ring bond stretch entries
// maxnel  maximum number of electronegativity bond corrections
// bcon    force constant parameters for harmonic bond stretch
// bcon5   force constant parameters for 5-ring bond stretch
// bcon4   force constant parameters for 4-ring bond stretch
// bcon3   force constant parameters for 3-ring bond stretch
// blen    bond length parameters for harmonic bond stretch
// blen5   bond length parameters for 5-ring bond stretch
// blen4   bond length parameters for 4-ring bond stretch
// blen3   bond length parameters for 3-ring bond stretch
// dlen    electronegativity bond length correction parameters
// kb      string of atom classes for harmonic bond stretch
// kb5     string of atom classes for 5-ring bond stretch
// kb4     string of atom classes for 4-ring bond stretch
// kb3     string of atom classes for 3-ring bond stretch
// kel     string of atom classes for electronegativity corrections

MDQC_EXTERN int maxnb;
MDQC_EXTERN int maxnb5;
MDQC_EXTERN int maxnb4;
MDQC_EXTERN int maxnb3;
MDQC_EXTERN int maxnel;
MDQC_EXTERN MDQCArray<real> bcon;
MDQC_EXTERN MDQCArray<real> bcon5;
MDQC_EXTERN MDQCArray<real> bcon4;
MDQC_EXTERN MDQCArray<real> bcon3;
MDQC_EXTERN MDQCArray<real> blen;
MDQC_EXTERN MDQCArray<real> blen5;
MDQC_EXTERN MDQCArray<real> blen4;
MDQC_EXTERN MDQCArray<real> blen3;
MDQC_EXTERN MDQCArray<real> dlen;
MDQC_EXTERN MDQCArray<std::string> kb;
MDQC_EXTERN MDQCArray<std::string> kb5;
MDQC_EXTERN MDQCArray<std::string> kb4;
MDQC_EXTERN MDQCArray<std::string> kb3;
MDQC_EXTERN MDQCArray<std::string> kel;
}
