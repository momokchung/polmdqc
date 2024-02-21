// Author: Moses KJ Chung
// Year:   2024

#include "tetrahedron.h"

namespace polmdqc
{
///////////////////////////////////////////////////////////////////
//                                                               //
//  tetrahedron  --  tetrahedron used in Delaunay/Alpha complex  //
//                                                               //
///////////////////////////////////////////////////////////////////

// "tetrahedron" class characterizes the tetrahedron used in Delaunay/Alpha
// complex theory; here we define the constructor, initializer, and destructor

Tetrahedron::Tetrahedron() {
    for (int i = 0; i < 4; i++) {
        this->vertices[i] = -1;
        this->neighbors[i] = -1;
        this->nindex[i] = -1;
    }
    std::bitset<8> b(std::string("00000000"));
    this->info = b;
}

void Tetrahedron::init() {
    for (int i = 0; i < 4; i++) {
        this->vertices[i] = -1;
        this->neighbors[i] = -1;
        this->nindex[i] = -1;
    }
    std::bitset<8> b(std::string("00000000"));
    this->info = b;
}

Tetrahedron::~Tetrahedron() {}
}
