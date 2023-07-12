////////////////////////////////////////////////////////
//                                                    //
//  mechanic.cpp  --  initialize molecular mechanics  //
//                                                    //
////////////////////////////////////////////////////////

// "mechanic" sets up needed parameters for the potential energy
// calculation and reads in many of the user selectable options


#include "active.h"
#include "angles.h"
#include "attach.h"
#include "bonds.h"
#include "inform.h"
#include "mechanic.h"

void mechanic()
{
    // set the bonded connectivity lists and active atoms
    attach();
    active();

    // find bonds, angles, torsions, bitorsions and small rings
    bonds();
    angles();
    // torsions();
    // bitors();
    // rings();
}
