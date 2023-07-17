////////////////////////////////////////////////////////
//                                                    //
//  mechanicqm.cpp  --  initialize quantum mechanics  //
//                                                    //
////////////////////////////////////////////////////////


#include "boys.h"
#include "katomsqm.h"
#include "kbasis.h"
#include "kgbs.h"
#include "kprim.h"
#include "kworker.h"
#include "mechanic.h"
#include <string>
#include <iostream>

namespace mechanicqm
{
//////////////////////////////////////////////////
//                                              //
//  mechanic  --  initialize quantum mechanics  //
//                                              //
//////////////////////////////////////////////////

void mechanic(std::string fileName)
{
    // read xyz file
    atoms::readxyz(fileName);

    // read gbs file
    gbs::readgbs(gbs::basisName);

    // set up basis
    basis::kbasis();

    // set up primitives
    prim::kprim();

    // read in boys function coefficients
    boys::initBoys();

    // allocate worker arrays
    worker::allocateWorker();
}


///////////////////////////////////////////////
//                                           //
//  cleanup  --  clean up quantum mechanics  //
//                                           //
///////////////////////////////////////////////

void cleanup()
{
    // cleanup boys function
    boys::cleanupBoys();
}
}
