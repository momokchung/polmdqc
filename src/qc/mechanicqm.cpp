// Author: Moses KJ Chung
// Year:   2024

#include "clusterqm.h"
#include "fatal.h"
#include "getbasis.h"
#include "inform.h"
#include "katomqm.h"
#include "mechanicqm.h"

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  mechanicqm  --  initialize quantum mechanics  //
//                                                //
////////////////////////////////////////////////////

// "mechanicqm" sets up needed parameters and basis for
// quantum mechanics calculation and reads in many of
// the user selectable options

void mechanicqm()
{
    // get the default basis set parameters
    getbasis();

    // assign atomic information
    katomqm();

    // set the atom groups
    clusterqm();

    // quit if essential parameter information is missing
    if (informAbort) {
        printf("\n MECHANICQM  --  Some Required Potential Energy Parameters are Undefined\n");
        fatal();
    }
}
}
