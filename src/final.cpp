////////////////////////////////////////////////////////
//                                                    //
//  final.cpp  --  final actions before program exit  //
//                                                    //
////////////////////////////////////////////////////////

// "final" performs any final program actions such as deallocation
// of global memory, prints a status message, and then pauses if
// necessary to avoid closing the execution window


#include "final.h"
#include "inform.h"

void final()
{
    // print a final status message before exiting PolQCMD
    if (debug) {
        printf("\n PolQCMD is Exiting following Normal Termination\n");
    }
}
