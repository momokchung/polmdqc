// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "precision.h"
#include <bitset>

namespace polmdqc
{
///////////////////////////////////////////////////////////////////
//                                                               //
//  tetrahedron  --  tetrahedron used in Delaunay/Alpha complex  //
//                                                               //
///////////////////////////////////////////////////////////////////

// vertices    the four vertices defining the tetrahedron
// neighbors   the four neighbors (-1 if no neighbour)
// info        0=orientation, 1=status, 2-5=surface info
// info_edge   
// nindex      vertex index of tetrahedron with common face
//             e.g. if tetra is (a,b,c,d), face (b,c,d) is also a
//             face of tetrahedron (b,c,d,e). (nindex=0,1,2,3)

class Tetrahedron {
public:
    int vertices[4];
    int neighbors[4];
    std::bitset<8> info;
    int info_edge[6];
    int nindex[4];

    Tetrahedron();
    ~Tetrahedron();

    void init();
};
}
