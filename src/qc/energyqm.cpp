// Author: Moses KJ Chung
// Year:   2023

#include "energyqm.h"
#include "hartree.h"

namespace energyqm
{
/////////////////////////////////////////
//                                     //
//  energy  --  compute system energy  //
//                                     //
/////////////////////////////////////////

void energy()
{
    // Run RHF for now. Will add other options later
    hartree::rhf();
}
}
