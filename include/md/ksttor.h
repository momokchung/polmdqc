// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  ksttor  --  stretch-torsion forcefield parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// maxnbt   maximum number of stretch-torsion parameter entries
// btcon    torsional amplitude parameters for stretch-torsion
// kbt      string of atom classes for stretch-torsion terms

MDQC_EXTERN int maxnbt;
MDQC_EXTERN MDQCArray2D<real,9> btcon;
MDQC_EXTERN MDQCArray<std::string> kbt;
}
