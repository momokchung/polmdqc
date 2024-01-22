// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

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
MDQC_EXTERN MDQCArray<real> kpep;
MDQC_EXTERN MDQCArray<real> prepep;
MDQC_EXTERN MDQCArray<real> dmppep;
MDQC_EXTERN MDQCArray3D<real,3,3> polscale;
MDQC_EXTERN MDQCArray3D<real,3,3> polinv;
MDQC_EXTERN MDQCArray<bool> lpep;
}
