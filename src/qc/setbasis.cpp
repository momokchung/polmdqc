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
    ngbs.allocate(maxele);
    namegbs.allocate(maxele);
    typgbs.allocate(maxele);
    amgbs.allocate(maxele);
    nprimgbs.allocate(maxele);
    scalegbs.allocate(maxele);
    coeffgbs.allocate(maxele);
    expgbs.allocate(maxele);
}
}
