///////////////////////////////////////////////////////////
//                                                       //
//  ktorsn.h  --  torsional angle forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

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


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnt;
QCMD_EXTERN int maxnt5;
QCMD_EXTERN int maxnt4;
QCMD_EXTERN std::vector<std::vector<double>> t1;
QCMD_EXTERN std::vector<std::vector<double>> t2;
QCMD_EXTERN std::vector<std::vector<double>> t3;
QCMD_EXTERN std::vector<std::vector<double>> t4;
QCMD_EXTERN std::vector<std::vector<double>> t5;
QCMD_EXTERN std::vector<std::vector<double>> t6;
QCMD_EXTERN std::vector<std::vector<double>> t15;
QCMD_EXTERN std::vector<std::vector<double>> t25;
QCMD_EXTERN std::vector<std::vector<double>> t35;
QCMD_EXTERN std::vector<std::vector<double>> t45;
QCMD_EXTERN std::vector<std::vector<double>> t55;
QCMD_EXTERN std::vector<std::vector<double>> t65;
QCMD_EXTERN std::vector<std::vector<double>> t14;
QCMD_EXTERN std::vector<std::vector<double>> t24;
QCMD_EXTERN std::vector<std::vector<double>> t34;
QCMD_EXTERN std::vector<std::vector<double>> t44;
QCMD_EXTERN std::vector<std::vector<double>> t54;
QCMD_EXTERN std::vector<std::vector<double>> t64;
QCMD_EXTERN std::vector<std::string> kt;
QCMD_EXTERN std::vector<std::string> kt5;
QCMD_EXTERN std::vector<std::string> kt4;
