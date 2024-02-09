// Author: Moses KJ Chung
// Year:   2024

#include "files.h"
#include "final.h"
#include "getcartqm.h"
#include "inform.h"
#include "initialqm.h"
#include "mechanicqm.h"
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
    getcartqm(ffile);
    mechanicqm();

    // perform any final tasks before program exit
    ffile.close();
    if (!test) final();
}
}
