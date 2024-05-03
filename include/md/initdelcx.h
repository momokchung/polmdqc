// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "alfatom.h"
#include "tetrahedron.h"
#include "vertex.h"
#include <queue>
#include <stack>
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////////////
//                                                           //
//  initdelcx  --  initialize/set up Delaunay triangulation  //
//                                                           //
///////////////////////////////////////////////////////////////

void initdelcx(int natoms, AlfAtom* alfatoms, std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra,
    std::queue<std::pair<int,int>>& link_facet, std::queue<std::pair<int,int>>& link_index, std::stack<int>& free, std::vector<int>& kill);
}
