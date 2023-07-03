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
#include <string>

const int maxele = 112;
extern double atmass[maxele];
extern double vdwrad[maxele];
extern double covrad[maxele];
extern std::string elemnt[maxele];
