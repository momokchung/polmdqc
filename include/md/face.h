// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "precision.h"
#include <bitset>

namespace polmdqc
{
////////////////////////////////////////////
//                                        //
//  face  --  face used in Alpha complex  //
//                                        //
////////////////////////////////////////////

class Face {
public:
    int vertices[3];
    int edges[3];
    real gamma;

    Face() {}

    Face(int i, int j, int k, int e1, int e2, int e3, real S);

    ~Face();
};
}
