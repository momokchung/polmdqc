// Author: Moses KJ Chung
// Year:   2023

#include "energyqm.h"
#include "init.h"
#include "mechanicqm.h"
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////////////
//                                                           //
//  tinkerqm  --  main function to run TinkerQM calculation  //
//                                                           //
///////////////////////////////////////////////////////////////

void tinkerqm(int argc, char** argv)
{
    // initialize
    init::init(argv);

    // take xyz file as first argument
    std::string fileName = argv[1];

    // set up quantum mechanics
    mechanicqm::mechanic(fileName);

    // add logic as to which function to run
    // for now we will just do HF
    energyqm::energy();

    // clean up quantum mechanics
    mechanicqm::cleanup();
}
}
