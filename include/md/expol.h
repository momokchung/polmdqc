// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  expol  --  exch-polarization in current structure  //
//                                                     //
/////////////////////////////////////////////////////////

// nexpol     total number of exch polarization sites in system
// kpep       exchange polarization spring constant at each site
// prepep     exchange polarization prefactor at each site
// dmppep     exchange polarization damping alpha at each site
// polscale   scale matrix for use in exchange polarization
// polinv     scale matrix inverse for exchange polarization
// lpep       flag to use exchange polarization at each site

MDQC_EXTERN int nexpol;
MDQC_EXTERN std::vector<real> kpep;
MDQC_EXTERN std::vector<real> prepep;
MDQC_EXTERN std::vector<real> dmppep;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> polscale;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> polinv;
MDQC_EXTERN std::vector<bool> lpep;
}
