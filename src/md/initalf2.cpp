// Author: Moses KJ Chung
// Year:   2024

#include "hilbert.h"
#include "initalf2.h"

namespace polmdqc
{
//////////////////////////////////////////
//                                      //
//  initalf2  --  initialize AlphaMol2  //
//                                      //
//////////////////////////////////////////

// "initalf2" initialize variables used in AlphaMol2

void initalf2()
{
    // initialize hilbert permutation table
    int ndim = 3;
    initHilbert(ndim);
}
}
