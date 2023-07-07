////////////////////////////////////////////////////
//                                                //
//  kinetic.h  --  compute kinetic energy matrix  //
//                                                //
////////////////////////////////////////////////////


#pragma once
#include "init.h"
#include <vector>

namespace kinetic
{
extern std::vector<std::vector<real>> cartKE;
extern std::vector<std::vector<real>> sphKE;

void kineticOS();
}
