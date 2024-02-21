// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "precision.h"
#include <bitset>
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  vertex  --  vertex used in Delaunay/Alpha complex  //
//                                                     //
/////////////////////////////////////////////////////////

// r        radius of vertex
// x        x Cartesian coordinate
// y        y Cartesian coordinate
// z        z Cartesian coordinate
// w        weight of vertex
// coefs    coefficient for surface
// coefv    coefficient for volume
// coefm    coefficient for mean curvature
// coefg    coefficient for gaussian curvature
// gamma    fraction of angle
// info     
// status   status of vertex (true if true vertex, false if bogus or infinite)

class Vertex {
    public:
        real r;
        real x,y,z;
        real w;
        real coefs, coefv, coefm, coefg;
        real gamma;

        std::bitset<8> info;
        bool status;

        Vertex() {}

        Vertex(real x, real y, real z, real r, real coefs, real coefv, real coefm, real coefg);

        ~Vertex();
};
}
