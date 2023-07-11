////////////////////////////////////////////////////////
//                                                    //
//  mechanic.cpp  --  initialize molecular mechanics  //
//                                                    //
////////////////////////////////////////////////////////

// "mechanic" sets up needed parameters for the potential energy
// calculation and reads in many of the user selectable options


#include "active.h"
#include "attach.h"
#include "inform.h"
#include "mechanic.h"

void mechanic()
{
    // set the bonded connectivity lists and active atoms
    attach();
    active();
}
