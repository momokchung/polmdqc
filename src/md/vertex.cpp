// Author: Moses KJ Chung
// Year:   2024

#include "vertex.h"

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  vertex  --  vertex used in Delaunay/Alpha complex  //
//                                                     //
/////////////////////////////////////////////////////////

// "vertex" class characterizes the vertices used in Delaunay/Alpha
// complex theory; here we define the constructor and destructor

Vertex::Vertex(real x, real y, real z, real r, real coefs, real coefv, real coefm, real coefg)
{
    this->coord[0] = x;
    this->coord[1] = y;
    this->coord[2] = z;
    this->r = r;
    this->coefs = coefs;
    this->coefv = coefv;
    this->coefm = coefm;
    this->coefg = coefg;

    std::bitset<8> b(std::string("00000000"));
    this->info = b;
    this->info[1] = 1;

    this->w = x*x + y*y + z*z - r*r;

    this->gamma = 0;
}

Vertex::~Vertex() {}
}
