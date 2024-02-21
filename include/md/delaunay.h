// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"
#include "tetrahedron.h"
#include "vertex.h"
#include <queue>
#include <stack>
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////
//                                              //
//  delaunay  --  Alpha shape theory variables  //
//                                              //
//////////////////////////////////////////////////

// vertices     vertices
// tetra        tetrahedrons
// link_facet   
// link_index   
// free         
// kill         

MDQC_EXTERN std::vector<Vertex> vertices;
MDQC_EXTERN std::vector<Tetrahedron> tetra;
MDQC_EXTERN std::queue<std::pair<int,int>> link_facet;
MDQC_EXTERN std::queue<std::pair<int,int>> link_index;
MDQC_EXTERN std::stack<int> free;
MDQC_EXTERN std::vector<int> kill;
}
