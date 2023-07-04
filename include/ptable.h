////////////////////////////////////////////////////////////
//                                                        //
//  ptable.h  --  symbols and info for chemical elements  //
//                                                        //
////////////////////////////////////////////////////////////

// maxele   maximum number of elements from periodic table
// atmass   standard atomic weight for each chemical element
// vdwrad   van der Waals radius for each chemical element
// covrad   covalent radius for each chemical element
// elemnt   atomic symbol for each chemical element


#pragma once
#include "macro.h"
#include <string>

const int maxele = 112;
QCMD_EXTERN double atmass[maxele];
QCMD_EXTERN double vdwrad[maxele];
QCMD_EXTERN double covrad[maxele];
QCMD_EXTERN std::string elemnt[maxele];
