// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  ktorsn  --  torsional angle forcefield parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// maxnt    maximum number of torsional angle parameter entries
// maxnt5   maximum number of 5-membered ring torsion entries
// maxnt4   maximum number of 4-membered ring torsion entries
// t1       torsional parameters for standard 1-fold rotation
// t2       torsional parameters for standard 2-fold rotation
// t3       torsional parameters for standard 3-fold rotation
// t4       torsional parameters for standard 4-fold rotation
// t5       torsional parameters for standard 5-fold rotation
// t6       torsional parameters for standard 6-fold rotation
// t15      torsional parameters for 1-fold rotation in 5-ring
// t25      torsional parameters for 2-fold rotation in 5-ring
// t35      torsional parameters for 3-fold rotation in 5-ring
// t45      torsional parameters for 4-fold rotation in 5-ring
// t55      torsional parameters for 5-fold rotation in 5-ring
// t65      torsional parameters for 6-fold rotation in 5-ring
// t14      torsional parameters for 1-fold rotation in 4-ring
// t24      torsional parameters for 2-fold rotation in 4-ring
// t34      torsional parameters for 3-fold rotation in 4-ring
// t44      torsional parameters for 4-fold rotation in 4-ring
// t54      torsional parameters for 5-fold rotation in 4-ring
// t64      torsional parameters for 6-fold rotation in 4-ring
// kt       string of atom classes for torsional angles
// kt5      string of atom classes for 5-ring torsions
// kt4      string of atom classes for 4-ring torsions

MDQC_EXTERN int maxnt;
MDQC_EXTERN int maxnt5;
MDQC_EXTERN int maxnt4;
MDQC_EXTERN std::vector<std::vector<real>> t1;
MDQC_EXTERN std::vector<std::vector<real>> t2;
MDQC_EXTERN std::vector<std::vector<real>> t3;
MDQC_EXTERN std::vector<std::vector<real>> t4;
MDQC_EXTERN std::vector<std::vector<real>> t5;
MDQC_EXTERN std::vector<std::vector<real>> t6;
MDQC_EXTERN std::vector<std::vector<real>> t15;
MDQC_EXTERN std::vector<std::vector<real>> t25;
MDQC_EXTERN std::vector<std::vector<real>> t35;
MDQC_EXTERN std::vector<std::vector<real>> t45;
MDQC_EXTERN std::vector<std::vector<real>> t55;
MDQC_EXTERN std::vector<std::vector<real>> t65;
MDQC_EXTERN std::vector<std::vector<real>> t14;
MDQC_EXTERN std::vector<std::vector<real>> t24;
MDQC_EXTERN std::vector<std::vector<real>> t34;
MDQC_EXTERN std::vector<std::vector<real>> t44;
MDQC_EXTERN std::vector<std::vector<real>> t54;
MDQC_EXTERN std::vector<std::vector<real>> t64;
MDQC_EXTERN std::vector<std::string> kt;
MDQC_EXTERN std::vector<std::string> kt5;
MDQC_EXTERN std::vector<std::string> kt4;
}
