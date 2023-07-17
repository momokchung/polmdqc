///////////////////////////////////////////////////////////
//                                                       //
//  kbonds.h  --  bond stretching forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// maxnb   maximum number of bond stretch parameter entries
// maxnb5  maximum number of 5-membered ring bond stretch entries
// maxnb4  maximum number of 4-membered ring bond stretch entries
// maxnb3  maximum number of 3-membered ring bond stretch entries
// maxnel  maximum number of electronegativity bond corrections
// bcon    force constant parameters for harmonic bond stretch
// bcon5   force constant parameters for 5-ring bond stretch
// bcon4   force constant parameters for 4-ring bond stretch
// bcon3   force constant parameters for 3-ring bond stretch
// blen    bond length parameters for harmonic bond stretch
// blen5   bond length parameters for 5-ring bond stretch
// blen4   bond length parameters for 4-ring bond stretch
// blen3   bond length parameters for 3-ring bond stretch
// dlen    electronegativity bond length correction parameters
// kb      string of atom classes for harmonic bond stretch
// kb5     string of atom classes for 5-ring bond stretch
// kb4     string of atom classes for 4-ring bond stretch
// kb3     string of atom classes for 3-ring bond stretch
// kel     string of atom classes for electronegativity corrections


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxnb;
QCMD_EXTERN int maxnb5;
QCMD_EXTERN int maxnb4;
QCMD_EXTERN int maxnb3;
QCMD_EXTERN int maxnel;
QCMD_EXTERN std::vector<double> bcon;
QCMD_EXTERN std::vector<double> bcon5;
QCMD_EXTERN std::vector<double> bcon4;
QCMD_EXTERN std::vector<double> bcon3;
QCMD_EXTERN std::vector<double> blen;
QCMD_EXTERN std::vector<double> blen5;
QCMD_EXTERN std::vector<double> blen4;
QCMD_EXTERN std::vector<double> blen3;
QCMD_EXTERN std::vector<double> dlen;
QCMD_EXTERN std::vector<std::string> kb;
QCMD_EXTERN std::vector<std::string> kb5;
QCMD_EXTERN std::vector<std::string> kb4;
QCMD_EXTERN std::vector<std::string> kb3;
QCMD_EXTERN std::vector<std::string> kel;
