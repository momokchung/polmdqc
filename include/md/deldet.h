// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "precision.h"

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  deldet  --  determinant for Delaunay triangulation  //
//                                                      //
//////////////////////////////////////////////////////////

void minor2(real a11, real a21, int& res, real eps=1e-10);
void minor3(real a11, real a12, real a21, real a22, real a31, real a32, int& res, real eps=1e-10);
void minor4(real* coord_a, real* coord_b, real* coord_c, real* coord_d, int& res, real eps=1e-10);
void minor5(real* coord_a, real r1, real* coord_b, real r2, real* coord_c,
    real r3, real* coord_d, real r4, real* coord_e, real r5, int& res, real eps=1e-10);
}
