////////////////////////////////////////////
//                                        //
//  eri.h  --  electron repulsion matrix  //
//                                        //
////////////////////////////////////////////


#pragma once
#include "init.h"
#include <vector>

namespace eri
{
extern std::vector<real> cartERI;
extern std::vector<real> sphERI;

void eriOS();
}
