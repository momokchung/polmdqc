// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "alfatom.h"

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  alphamol  --  compute surface area and volume  //
//                                                 //
/////////////////////////////////////////////////////

void alphamol(int natoms, AlfAtom* alfatoms, real& wsurf, real& wvol, real& wmean, real& wgauss,
    real* surf, real* vol, real* mean, real* gauss,
    real* dsurf, real* dvol, real* dmean, real* dgauss, bool deriv);
}
