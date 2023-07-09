///////////////////////////////////////////////////////
//                                                   //
//  fatal.cpp  --  terminate the program abnormally  //
//                                                   //
///////////////////////////////////////////////////////

// "fatal" terminates execution due to a user request, a severe
// error or some other nonstandard condition


#include "fatal.h"
#include "final.h"
#include <cstdlib>
#include <stdio.h>

void fatal()
{
    // print a final warning message, then do final cleanup
    printf("\n\n PolQCMD is Unable to Continue; Terminating the Current Calculation\n");
    final();
    exit(1);
}
