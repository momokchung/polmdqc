// Author: Moses KJ Chung
// Year:   2023

#include "fatal.h"
#include "final.h"
#include <cstdlib>
#include <stdio.h>

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  fatal  --  terminate the program abnormally  //
//                                               //
///////////////////////////////////////////////////

// "fatal" terminates execution due to a user request, a severe
// error or some other nonstandard condition

void fatal()
{
    // print a final warning message, then do final cleanup
    printf("\n PolMDQC is Unable to Continue; Terminating the Current Calculation\n");
    final();
    exit(1);
}
}
