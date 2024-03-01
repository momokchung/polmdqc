// Author: Moses KJ Chung
// Year:   2024

#include "face.h"

namespace polmdqc
{
////////////////////////////////////////////
//                                        //
//  face  --  face used in Alpha complex  //
//                                        //
////////////////////////////////////////////

// "face" class characterizes the face used in Alpha
// complex theory

Face::Face(int i, int j, int k, int e1, int e2, int e3, real S) {
    vertices[0] = i;
    vertices[1] = j;
    vertices[2] = k;
    edges[0] = e1;
    edges[1] = e2;
    edges[2] = e3;
    gamma = S;
}

Face::~Face() {}
}
