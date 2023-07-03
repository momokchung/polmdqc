///////////////////////////////////////////////////////////////////
//                                                               //
//  tinkerqm.cpp  --  main function to run TinkerQM calculation  //
//                                                               //
///////////////////////////////////////////////////////////////////


#include "energy.h"
#include "mechanic.h"
#include <string>

int main(int argc, char** argv)
{
    // take xyz file as first argument
    std::string fileName = argv[1];

    // set up quantum mechanics
    mechanic::mechanic(fileName);

    // add logic as to which function to run
    // for now we will just do HF
    energy::energy();

    // clean up quantum mechanics
    mechanic::cleanup();
}
