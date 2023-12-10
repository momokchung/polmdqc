// Author: Moses KJ Chung
// Year:   2023

#include "final.h"
#include "inform.h"

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  final  --  final actions before program exit  //
//                                                //
////////////////////////////////////////////////////

// "final" performs any final program actions such as deallocation
// of global memory, prints a status message, and then pauses if
// necessary to avoid closing the execution window

void final()
{
    // print a final status message before exiting PolMDQC
    if (debug) {
        printf("\n PolMDQC is Exiting following Normal Termination\n");
    }
}
}