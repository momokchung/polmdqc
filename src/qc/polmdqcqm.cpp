// Author: Moses KJ Chung
// Year:   2024

#include "initialqm.h"
#include "polmdqcqm.h"

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  polmdqcqm  --  energy partitioning and analysis  //
//                                                   //
///////////////////////////////////////////////////////

// "polmdqcqm" computes the quantum mechanical energy

void polmdqcqm(int argc, char** argv)
{
    // set up the structure and mechanics calculation
    initialqm(argc, argv);
    // getcartqm(ffile);
    // mechanicqm();
}
}
