// Author: Moses KJ Chung
// Year:   2024

#include "edge.h"

namespace polmdqc
{
////////////////////////////////////////////
//                                        //
//  edge  --  edge used in Alpha complex  //
//                                        //
////////////////////////////////////////////

// "edge" class characterizes the edge used in Alpha
// complex theory

Edge::Edge(int i, int j) {
    vertices[0] = i;
    vertices[1] = j;
    gamma = 0.0;
    sigma = 0.0;
}

Edge::~Edge() {}
}
