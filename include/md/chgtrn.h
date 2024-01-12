// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  chgtrn  --  charge transfer in current structure  //
//                                                    //
////////////////////////////////////////////////////////

// nct       total number of dispersion sites in the system
// chgct     charge for charge transfer at each multipole site
// dmpct     charge transfer damping factor at each multipole site

MDQC_EXTERN int nct;
MDQC_EXTERN std::vector<real> chgct;
MDQC_EXTERN std::vector<real> dmpct;
}
