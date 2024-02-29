// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "tetrahedron.h"

namespace polmdqc
{
////////////////////////////////////////
//                                    //
//  findedge  --  find index of edge  //
//                                    //
////////////////////////////////////////

// "findedge" finds the index of the edge
// given two vertices of a tetrahedron

inline int findedge(Tetrahedron t, int i1, int j1)
{
    int ipair;

    if (i1==t.vertices[0]) {
        if (j1==t.vertices[1]) ipair = 5;
        else if (j1==t.vertices[2]) ipair = 4;
        else ipair = 3;
    }
    else if (i1==t.vertices[1]) {
        if (j1==t.vertices[2]) ipair = 2;
        else ipair = 1;
    }
    else {
        ipair = 0;
    }
    return ipair;
}
}
