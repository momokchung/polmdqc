// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "kmulti.h"
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  mpole  --  atomic multipoles in current structure  //
//                                                     //
/////////////////////////////////////////////////////////

// maxpole   max components (monopole=1,dipole=4,quadrupole=13)
// npole     total number of multipole sites in the system
// ipole     number of the atom for each multipole site
// polsiz    number of multipole components for each atom
// pollist   multipole site for each atom (-1=no multipole)
// zaxis     number of the z-axis defining atom for each atom (starting index=1)
// xaxis     number of the x-axis defining atom for each atom (starting index=1)
// yaxis     number of the y-axis defining atom for each atom (starting index=1)
// mono0     original atomic monopole values for charge flux
// mscale    multipole exclusion coefficient scale
// pole      local frame Cartesian multipoles for each atom
// rpole     global frame Cartesian multipoles for each atom
// tem       multipole torque for each atom
// polaxe    local coordinate frame type for each atom

constexpr int maxpole = 13;
MDQC_EXTERN int npole;
MDQC_EXTERN MDQCArray<int> ipole;
MDQC_EXTERN MDQCArray<int> polsiz;
MDQC_EXTERN MDQCArray<int> pollist;
MDQC_EXTERN MDQCArray<int> zaxis;
MDQC_EXTERN MDQCArray<int> xaxis;
MDQC_EXTERN MDQCArray<int> yaxis;
MDQC_EXTERN MDQCArray<real> mono0;
MDQC_EXTERN MDQCArray<real> mscale;
MDQC_EXTERN MDQCArray2D<real,maxpole> pole;
MDQC_EXTERN MDQCArray2D<real,maxpole> rpole;
MDQC_EXTERN MDQCArray2D<real,3> tem;
MDQC_EXTERN MDQCArray<LocalFrame> polaxe;
}
