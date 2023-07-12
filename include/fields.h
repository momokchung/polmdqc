//////////////////////////////////////////////////////////
//                                                      //
//  fields.h  --  molecular mechanics force field type  //
//                                                      //
//////////////////////////////////////////////////////////

// biotyp       force field atom type of each biopolymer type
// forcefield   string used to describe the current forcefield


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN std::vector<int> biotyp;
QCMD_EXTERN std::string forcefield;
