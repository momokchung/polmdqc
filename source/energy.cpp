/////////////////////////////////////////////
//                                         //
//  energy.cpp  --  compute system energy  //
//                                         //
/////////////////////////////////////////////


#include "energy.h"
#include "hartree.h"

namespace energy
{
//////////////////////////////////////////////
//                                          //
//  void energy  --  compute system energy  //
//                                          //
//////////////////////////////////////////////

void energy()
{
    // Run RHF for now. Will add other options later
    hartree::rhf();
}
}