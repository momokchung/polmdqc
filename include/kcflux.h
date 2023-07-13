//////////////////////////////////////////////////////////
//                                                      //
//  kcflux.h -- charge flux term forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// maxncfb   maximum number of bond stretch charge flux entries
// maxncfa   maximum number of angle bend charge flux entries
// cflb      charge flux over stretching of a bond length
// cfla      charge flux over bending of a bond angle
// cflab     charge flux over asymmetric bond within an angle
// kcfb      string of atom classes for bond stretch charge flux
// kcfa      string of atom classes for angle bend charge flux


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxncfb;
QCMD_EXTERN int maxncfa;
QCMD_EXTERN std::vector<double> cflb;
QCMD_EXTERN std::vector<std::vector<double>> cfla;
QCMD_EXTERN std::vector<std::vector<double>> cflab;
QCMD_EXTERN std::vector<std::string>> kcfb;
QCMD_EXTERN std::vector<std::string>> kcfa;
