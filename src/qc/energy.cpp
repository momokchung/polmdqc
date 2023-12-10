// Author: Moses KJ Chung
// Year:   2023

#include "energy.h"
#include "hartree.h"

namespace energy
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
