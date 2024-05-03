// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "tetrahedron.h"
#include "vertex.h"
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////
//                                          //
//  alfcx  --  build alpha complex from DT  //
//                                          //
//////////////////////////////////////////////

void alfcx(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra, real alpha);
}
