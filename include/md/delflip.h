// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "tetrahedron.h"
#include "vertex.h"
#include <queue>
#include <stack>
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  delflip  --  flips used in Delaunay triangulation  //
//                                                     //
/////////////////////////////////////////////////////////

void flip_1_4(std::vector<Tetrahedron>& tetra, int ipoint, int itetra, int& tetra_last,
    std::queue<std::pair<int,int>>& link_facet, std::queue<std::pair<int,int>>& link_index, std::stack<int>& free, std::vector<int>& kill);
void flip(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra,
    std::queue<std::pair<int,int>>& link_facet, std::queue<std::pair<int,int>>& link_index, std::stack<int>& free, std::vector<int>& kill);
}
