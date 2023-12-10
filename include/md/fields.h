// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  fields  --  molecular mechanics force field type  //
//                                                    //
////////////////////////////////////////////////////////

// biotyp       force field atom type of each biopolymer type
// forcefield   string used to describe the current forcefield

MDQC_EXTERN std::vector<int> biotyp;
MDQC_EXTERN std::string forcefield;
}
