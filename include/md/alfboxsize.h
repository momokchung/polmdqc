// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "alfatom.h"

namespace polmdqc
{
/////////////////////////////////////////////////////////////
//                                                         //
//  alfboxsize  --  find smallest and largest coordinates  //
//                                                         //
/////////////////////////////////////////////////////////////

void alfboxsize(AlfAtom* alfatoms, int size, real& xmin, real& ymin, real& zmin, real& xmax, real& ymax, real& zmax, real& rmax);
}
