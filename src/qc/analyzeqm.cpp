// Author: Moses KJ Chung
// Year:   2024

#include "analyzeqm.h"
#include "files.h"
#include "final.h"
#include "getcartqm.h"
#include "inform.h"
#include "initialqm.h"
#include "mechanicqm.h"

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  analyzeqm  --  energy partitioning and analysis  //
//                                                   //
///////////////////////////////////////////////////////

// "analyzeqm" computes the quantum mechanical energy

void analyzeqm(int argc, char** argv)
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
