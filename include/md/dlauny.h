// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"
#include "edge.h"
#include "face.h"
#include "tetrahedron.h"
#include "vertex.h"
#include <queue>
#include <stack>
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  dlauny  --  Delaunay triangulation variables  //
//                                                //
////////////////////////////////////////////////////

// deleps       epsilon value
// delepsvol    epsilon volume value
// vertices     vertices
// tetra        tetrahedrons
// edges        edges
// faces        faces
// link_facet   
// link_index   
// free         
// kill         

constexpr real deleps = 1e-3;
constexpr real delepsvol = 1e-13;
MDQC_EXTERN std::vector<Vertex> vertices;
MDQC_EXTERN std::vector<Tetrahedron> tetra;
MDQC_EXTERN std::vector<Edge> edges;
MDQC_EXTERN std::vector<Face> faces;
MDQC_EXTERN std::queue<std::pair<int,int>> link_facet;
MDQC_EXTERN std::queue<std::pair<int,int>> link_index;
MDQC_EXTERN std::stack<int> free;
MDQC_EXTERN std::vector<int> kill;

constexpr int inf4_1[4] = {1, 1, 0, 0};
constexpr int sign4_1[4] = {-1, 1, 1, -1};
constexpr int inf4_2[4][4] = {
    { -1, 1, 2, 2},
    { 1, -1, 2, 2},
    { 2, 2, -1, 0},
    { 2, 2, 0, -1}
};
constexpr int sign4_2[4][4] = {
    { 0, 1, -1, 1},
    { -1, 0, 1, -1},
    { 1, -1, 0, 1},
    { -1, 1, -1, 0}
};
constexpr int sign4_3[4] = {-1, 1, -1, 1};
constexpr int inf5_2[4][4] = {
    { -1, 1, 0, 0},
    { 1, -1, 0, 0},
    { 0, 0, -1, 0},
    { 0, 0, 0, -1}
};
constexpr int sign5_2[4][4] = {
    { 0, -1, -1, 1},
    { 1, 0, -1, 1},
    { 1, 1, 0, 1},
    { -1, -1, -1, 0}
};
constexpr int inf5_3[4] = {0, 0, 2, 2};
constexpr int sign5_3[4] = {1, 1, -1, 1};
constexpr int order1[4][3] = {
    { 2, 1, 3},
    { 0, 2, 3},
    { 1, 0, 3},
    { 0, 1, 2}
};

constexpr int ord_rc[3][3] = {
    {0, 1, 2},
    {2, 0, 1},
    {1, 2, 0},
};

constexpr int order2[6][2] = {
    { 2, 3},
    { 3, 1},
    { 1, 2},
    { 0, 3},
    { 2, 0},
    { 0, 1}
};
constexpr int order3[6][2] = {
    { 0, 1},
    { 0, 2},
    { 0, 3},
    { 1, 2},
    { 1, 3},
    { 2, 3}
};
constexpr int idxList[4][3] = {
    { 0, 0, 0},
    { 0, 1, 1},
    { 1, 1, 2},
    { 2, 2, 2}
};
constexpr int table32[3][3] = {
    { 0, 1, 2},
    { 0, 2, 1},
    { 2, 0, 1}
};
constexpr int table32_2[3][2] = {
    { 0, 1},
    { 0, 2},
    { 1, 2}
};
constexpr int table41[3][3] = {
    { 1, 0, 2},
    { 0, 1, 2},
    { 0, 2, 1}
};
constexpr int table41_2[3][2] = {
    { 0, 0},
    { 1, 0},
    { 1, 1}
};
constexpr int order[3][2] = {
    { 1, 2},
    { 2, 0},
    { 0, 1}
};
constexpr int other[4][3] = {
    { 1, 2, 3},
    { 0, 2, 3},
    { 0, 1, 3},
    { 0, 1, 2}
};
constexpr int other2[4][4][2] = {
    {
        { -1, -1},
        { 2, 3},
        { 1, 3},
        { 1, 2}
    },
    {
        { 2, 3},
        { -1, -1},
        { 0, 3},
        { 0, 2}
    },
    {
        { 1, 3},
        { 0, 3},
        { -1, -1},
        { 0, 1}
    },
    {
        { 1, 2},
        { 0, 2},
        { 0, 1},
        { -1, -1}
    },
};
}
