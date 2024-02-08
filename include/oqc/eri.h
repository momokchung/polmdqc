// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "init.h"
#include <vector>

namespace eri
{
//////////////////////////////////////////
//                                      //
//  eri  --  electron repulsion matrix  //
//                                      //
//////////////////////////////////////////

extern std::vector<real> cartERI;
extern std::vector<real> sphERI;

void eriOS();
}
