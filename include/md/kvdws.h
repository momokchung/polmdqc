// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  kvdws  --  van der Waals term forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// rad      van der Waals radius parameter for each atom class
// eps      van der Waals well depth parameter for each atom class
// rad4     van der Waals radius parameter in 1-4 interactions
// eps4     van der Waals well depth parameter in 1-4 interactions
// reduct   van der Waals reduction factor for each atom class

MDQC_EXTERN std::vector<double> rad;
MDQC_EXTERN std::vector<double> eps;
MDQC_EXTERN std::vector<double> rad4;
MDQC_EXTERN std::vector<double> eps4;
MDQC_EXTERN std::vector<double> reduct;
}
