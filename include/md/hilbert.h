// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "dlauny2.h"

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  hilbert  --  sort points using a Hilbert Curve  //
//                                                  //
//////////////////////////////////////////////////////

void initHilbert(int ndim);

void sort3DHilbert(AlfAtom *alfatoms, int size, int e, int d, real xmin, real xmax, real ymin, real ymax, real zmin, real zmax, int depth);

void brioHilbert(AlfAtom *alfatoms, int size, real xmin, real xmax, real ymin, real ymax, real zmin, real zmax, int depth);
}
