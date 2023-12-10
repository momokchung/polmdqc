// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  titles  --  title for current molecular system  //
//                                                  //
//////////////////////////////////////////////////////

// ltitle   length in characters of the nonblank title string
// title    title used to describe the current structure

MDQC_EXTERN int ltitle;
MDQC_EXTERN std::string title;
}
