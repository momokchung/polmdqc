// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  chunks  --  PME grid spatial decomposition values  //
//                                                     //
/////////////////////////////////////////////////////////

// nchunk     total number of spatial regions for PME grid
// nchk1      number of spatial regions along the a-axis
// nchk2      number of spatial regions along the b-axis
// nchk3      number of spatial regions along the c-axis
// ngrd1      number of grid points per region along a-axis
// ngrd2      number of grid points per region along b-axis
// ngrd3      number of grid points per region along c-axis
// nlpts      PME grid points to the left of center point
// nrpts      PME grid points to the right of center point
// grdoff     offset for index into B-spline coefficients
// pmetable   PME grid spatial regions involved for each site

MDQC_EXTERN int nchunk;
MDQC_EXTERN int nchk1,nchk2,nchk3;
MDQC_EXTERN int ngrd1,ngrd2,ngrd3;
MDQC_EXTERN int nlpts,nrpts,grdoff;
MDQC_EXTERN std::vector<std::vector<int>> pmetable;
}
