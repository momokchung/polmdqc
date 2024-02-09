// Author: Moses KJ Chung
// Year:   2023

#include "darray.h"
#include "finalqm.h"
#include "mod.h"

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  finalqm  --  final actions before program exit  //
//                                                  //
//////////////////////////////////////////////////////

// "finalqm" performs any final program actions such as deallocation
// of global memory, prints a status message, and then pauses if
// necessary to avoid closing the execution window

void finalqm()
{
    // print a final status message before exiting PolMDQC
    if (debug) {
        printf("\n PolMDQC is Exiting following Normal Termination\n");
    }

    // deallocation of global arrays from module kgbs
    nprimgbs.deallocate();
    ngbs.deallocate();
    scalegbs.deallocate();
    namegbs.deallocate();
    amgbs.deallocate();
    typgbs.deallocate();
    coeffgbs.deallocate();
    expgbs.deallocate();
}
}
