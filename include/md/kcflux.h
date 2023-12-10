// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  kcflux -- charge flux term forcefield parameters  //
//                                                    //
////////////////////////////////////////////////////////

// maxncfb   maximum number of bond stretch charge flux entries
// maxncfa   maximum number of angle bend charge flux entries
// cflb      charge flux over stretching of a bond length
// cfla      charge flux over bending of a bond angle
// cflab     charge flux over asymmetric bond within an angle
// kcfb      string of atom classes for bond stretch charge flux
// kcfa      string of atom classes for angle bend charge flux

MDQC_EXTERN int maxncfb;
MDQC_EXTERN int maxncfa;
MDQC_EXTERN std::vector<double> cflb;
MDQC_EXTERN std::vector<std::vector<double>> cfla;
MDQC_EXTERN std::vector<std::vector<double>> cflab;
MDQC_EXTERN std::vector<std::string> kcfb;
MDQC_EXTERN std::vector<std::string> kcfa;
}
