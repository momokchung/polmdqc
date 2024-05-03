// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "face.h"
#include "tetrahedron.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  alfcxfaces  --  generates the list of faces  //
//                                               //
///////////////////////////////////////////////////

void alfcxfaces(std::vector<Tetrahedron>& tetra, std::vector<Face>& faces);
}
