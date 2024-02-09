// Author: Moses KJ Chung
// Year:   2024

#include "kgbs.h"
#include "ptable.h"
#include "setbasis.h"

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  setbasis  --  allocate basis set parameters  //
//                                               //
///////////////////////////////////////////////////

// "setbasis" allocates memory space for basis set parameters

void setbasis()
{
    // allocate basis set parameters 
    nprimgbs.allocate(maxele);
    ngbs.allocate(maxele);
    scalegbs.allocate(maxele);
    namegbs.allocate(maxele);
    amgbs.allocate(maxele);
    typgbs.allocate(maxele);
    coeffgbs.allocate(maxele);
    expgbs.allocate(maxele);
}
}
