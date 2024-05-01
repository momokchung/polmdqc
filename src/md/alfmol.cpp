// Author: Moses KJ Chung
// Year:   2024

#include "alfmol.h"
#include "alfp.h"
#include "alphamol.h"
#include "alphamol2.h"

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
    alphamol(deriv);
    // if (alfmeth == AlfMethod::AlphaMol) alphamol(deriv);
    // else if (alfmeth == AlfMethod::AlphaMol2) alphamol2(deriv);
}
}
