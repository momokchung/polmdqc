// Author: Moses KJ Chung
// Year:   2024

#include "alfmol.h"
#include "alfp.h"
#include "alphamol1.h"
#include "alphamol2.h"
#include "inform.h"
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
    clock_t start_s,stop_s;
    if (verbose) start_s = clock();
    initalfatm(deriv);
    if (verbose) {
        stop_s = clock();
        printf("\n Initalfatm compute time   : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
    }

    // run AlphaMol
    if (alfmeth == AlfMethod::AlphaMol) alphamol1(deriv);
    else if (alfmeth == AlfMethod::AlphaMol2) alphamol2(deriv);
}
}
