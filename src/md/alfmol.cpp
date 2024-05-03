// Author: Moses KJ Chung
// Year:   2024

#include "alfmol.h"
#include "alfp.h"
#include "alphamol.h"
#include "alphamol2.h"
#include "alphmol.h"
#include "dlauny2.h"
#include "initalfatm.h"

namespace polmdqc
{
/////////////////////////////////////////////
//                                         //
//  alfmol  --  run AlphaMol or AlphaMol2  //
//                                         //
/////////////////////////////////////////////

// "alfmol" runs AlphaMol or AlphaMol2 based on user input

void alfmol(bool deriv)
{
    // initialize alfatoms
    initalfatm();

    // run AlphaMol
    alphamol(alfatoms.size(), &(alfatoms[0]), wsurf, wvol, wmean, wgauss,
        surf.ptr(), vol.ptr(), mean.ptr(), gauss.ptr(),
        dsurf.ptr(), dvol.ptr(), dmean.ptr(), dgauss.ptr(), deriv);
    // if (alfmeth == AlfMethod::AlphaMol) alphamol(deriv);
    // else if (alfmeth == AlfMethod::AlphaMol2) alphamol2(deriv);
}
}
