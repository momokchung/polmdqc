// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "darray.h"
#include "precision.h"

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  alterchg  --  modification of partial charges  //
//                                                 //
/////////////////////////////////////////////////////

void alterchg();

void bndchg(real* pdelta);

void angchg(real* pdelta);
}
