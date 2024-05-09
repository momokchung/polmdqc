// Author: Moses KJ Chung
// Year:   2024

#include "alfp.h"
#include "vertex.h"
#include <cmath>

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
    real x1 = truncate_real(x, alfdigit);
    real y1 = truncate_real(y, alfdigit);
    real z1 = truncate_real(z, alfdigit);
    real r1 = truncate_real(r, alfdigit);
    this->coord[0] = x1;
    this->coord[1] = y1;
    this->coord[2] = z1;
    this->r = r1;
    this->coefs = coefs;
    this->coefv = coefv;
    this->coefm = coefm;
    this->coefg = coefg;

    std::bitset<8> b(std::string("00000000"));
    this->info = b;
    this->info[1] = 1;

    real t1 = (real) pow(10, alfdigit/2);
    real t2 = (real) pow(10, alfdigit);
    long long ival1= std::round(r1*t1);
    long long ival2 = -ival1*ival1;
    ival1 = std::round(x1*t1);
    ival2 += ival1*ival1;
    ival1 = std::round(y1*t1);
    ival2 += ival1*ival1;
    ival1 = std::round(z1*t1);
    ival2 += ival1*ival1;
    this->w = (real) ival2/t2;

    this->gamma = 0;
}

Vertex::~Vertex() {}

real Vertex::truncate_real(real x, int ndigit)
{
    real x_out,y;
    real fact;

    int mantissa;
    int digit;

    mantissa = (int) x;
    y = x - mantissa;

    x_out = mantissa;
    fact = 1;
    for (int i = 0; i < ndigit; i++) {
        fact *= 10;
        digit = (int) std::round(y*10);
        y = 10*(y-digit/10.0);
        x_out += digit/fact;
    }

    return x_out;
}
}
