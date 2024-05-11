// Author: Moses KJ Chung
// Year:   2024

#include "random.h"
#include <cmath>
#include <random>

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  ranvec  --  unit vector in random direction  //
//                                               //
///////////////////////////////////////////////////

// "ranvec" generates a unit vector in 3-dimensional
// space with uniformly distributed random orientation;
// zero is not included in the distribution

void ranvec(real vec[3])
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);
    std::uniform_real_distribution<real> distribution1(0.1, 1.0);
    std::uniform_real_distribution<real> distribution2(-1.0, -0.1);

    real r[3];
    for (int i = 0; i < 3; i++) {
        if (dis(gen) == 0) r[i] = distribution1(gen);
        else r[i] = distribution2(gen);
    }

    real mag = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

    // Normalize the vector
    vec[0] = r[0] / mag;
    vec[1] = r[1] / mag;
    vec[2] = r[2] / mag;
}

//////////////////////////////////////////////////////
//                                                  //
//  wiggle  --  random perturbation of coordinates  //
//                                                  //
//////////////////////////////////////////////////////

// "wiggle" applies a small random perturbation of coordinates
// to avoid numerical instability in geometric calculations for
// linear, planar and other symmetric structures

void wiggle(int n, real* x, real* y, real* z, real eps)
{
    real vec[3];
    for (int i = 0; i < n; i++) {
        ranvec(vec);
        x[i] += eps*vec[0];
        y[i] += eps*vec[1];
        z[i] += eps*vec[2];
    }
}
}
