////////////////////////////////////////////////////////////
//                                                        //
//  kitors.h  --  improper torsion forcefield parameters  //
//                                                        //
////////////////////////////////////////////////////////////

// maxnti   maximum number of improper torsion parameter entries
// ti1      torsional parameters for improper 1-fold rotation
// ti2      torsional parameters for improper 2-fold rotation
// ti3      torsional parameters for improper 3-fold rotation
// kti      string of atom classes for improper torsional parameters


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnti;
QCMD_EXTERN std::vector<std::vector<double>> ti1;
QCMD_EXTERN std::vector<std::vector<double>> ti2;
QCMD_EXTERN std::vector<std::vector<double>> ti3;
QCMD_EXTERN std::vector<std::string> kti;
