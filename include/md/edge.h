// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "precision.h"

namespace polmdqc
{
////////////////////////////////////////////
//                                        //
//  edge  --  edge used in Alpha complex  //
//                                        //
////////////////////////////////////////////

class Edge {
public:
    int vertices[2];
    real gamma;
    real sigma;
    real len,surf,vol; 
    real coefm1,coefm2,coefg1,coefg2;
    real dsurf,dvol,dmean,dgauss;

    Edge() {}

    Edge(int i, int j);

    ~Edge();
};
}
