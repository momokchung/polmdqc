//////////////////////////////////////////////////////////
//                                                      //
//  kangs.h  --  bond angle bend forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// maxna    maximum number of harmonic angle bend parameter entries
// maxna5   maximum number of 5-membered ring angle bend entries
// maxna4   maximum number of 4-membered ring angle bend entries
// maxna3   maximum number of 3-membered ring angle bend entries
// maxnap   maximum number of in-plane angle bend parameter entries
// maxnaf   maximum number of Fourier angle bend parameter entries
// acon     force constant parameters for harmonic angle bends
// acon5    force constant parameters for 5-ring angle bends
// acon4    force constant parameters for 4-ring angle bends
// acon3    force constant parameters for 3-ring angle bends
// aconp    force constant parameters for in-plane angle bends
// aconf    force constant parameters for Fourier angle bends
// ang      bond angle parameters for harmonic angle bends
// ang5     bond angle parameters for 5-ring angle bends
// ang4     bond angle parameters for 4-ring angle bends
// ang3     bond angle parameters for 3-ring angle bends
// angp     bond angle parameters for in-plane angle bends
// angf     phase shift angle and periodicity for Fourier bends
// ka       string of atom classes for harmonic angle bends
// ka5      string of atom classes for 5-ring angle bends
// ka4      string of atom classes for 4-ring angle bends
// ka3      string of atom classes for 3-ring angle bends
// kap      string of atom classes for in-plane angle bends
// kaf      string of atom classes for Fourier angle bends


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxna;
QCMD_EXTERN int maxna5;
QCMD_EXTERN int maxna4;
QCMD_EXTERN int maxna3;
QCMD_EXTERN int maxnap;
QCMD_EXTERN int maxnaf;
QCMD_EXTERN std::vector<double> acon;
QCMD_EXTERN std::vector<double> acon5;
QCMD_EXTERN std::vector<double> acon4;
QCMD_EXTERN std::vector<double> acon3;
QCMD_EXTERN std::vector<double> aconp;
QCMD_EXTERN std::vector<double> aconf;
QCMD_EXTERN std::vector<std::vector<double>> ang;
QCMD_EXTERN std::vector<std::vector<double>> ang5;
QCMD_EXTERN std::vector<std::vector<double>> ang4;
QCMD_EXTERN std::vector<std::vector<double>> ang3;
QCMD_EXTERN std::vector<std::vector<double>> angp;
QCMD_EXTERN std::vector<std::vector<double>> angf;
QCMD_EXTERN std::vector<std::string> ka;
QCMD_EXTERN std::vector<std::string> ka5;
QCMD_EXTERN std::vector<std::string> ka4;
QCMD_EXTERN std::vector<std::string> ka3;
QCMD_EXTERN std::vector<std::string> kap;
QCMD_EXTERN std::vector<std::string> kaf;
