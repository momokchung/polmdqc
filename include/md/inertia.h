// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "precision.h"

namespace polmdqc
{
/////////////////////////////////////////////////
//                                             //
//  inertia  --  principal moments of inertia  //
//                                             //
/////////////////////////////////////////////////

void inertia(int mode, int n, real* mass, real* x, real* y, real* z);
}
