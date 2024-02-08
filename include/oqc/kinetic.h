// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "init.h"
#include <vector>

namespace kinetic
{
//////////////////////////////////////////////////
//                                              //
//  kinetic  --  compute kinetic energy matrix  //
//                                              //
//////////////////////////////////////////////////

extern std::vector<std::vector<real>> cartKE;
extern std::vector<std::vector<real>> sphKE;

void kineticOS();
}
