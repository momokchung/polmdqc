// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "edge.h"
#include "tetrahedron.h"
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////
//                                          //
//  alfcxedges  --  generate list of edges  //
//                                          //
//////////////////////////////////////////////

void alfcxedges(std::vector<Tetrahedron>& tetra, std::vector<Edge>& edges);
}
