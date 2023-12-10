// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  urypot  --  Urey-Bradley functional form details  //
//                                                    //
////////////////////////////////////////////////////////

// cury       cubic coefficient in Urey-Bradley potential
// qury       quartic coefficient in Urey-Bradley potential
// ureyunit   convert Urey-Bradley energy to kcal/mole

MDQC_EXTERN double cury;
MDQC_EXTERN double qury;
MDQC_EXTERN double ureyunit;
}
