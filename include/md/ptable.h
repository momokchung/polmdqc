// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <map>
#include <string>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  ptable  --  symbols and info for chemical elements  //
//                                                      //
//////////////////////////////////////////////////////////

// maxele      maximum number of elements from periodic table
// atmass      standard atomic weight for each chemical element
// vdwrad      van der Waals radius for each chemical element
// covrad      covalent radius for each chemical element
// elemnt      atomic symbol for each chemical element
// symtoatmn   map from atomic symbol to atomic number

constexpr int maxele = 112;
MDQC_EXTERN real atmass[maxele];
MDQC_EXTERN real vdwrad[maxele];
MDQC_EXTERN real covrad[maxele];
MDQC_EXTERN std::string elemnt[maxele];
MDQC_EXTERN std::map<std::string, int> symtoatmn;
}
