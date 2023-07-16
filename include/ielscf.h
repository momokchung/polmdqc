/////////////////////////////////////////////////////////
//                                                     //
//  ielscf.h  --  extended Lagrangian induced dipoles  //
//                                                     //
/////////////////////////////////////////////////////////

// nfree_aux    total degrees of freedom for auxiliary dipoles
// tautemp_aux  time constant for auliliary Berendsen thermostat
// kelvin_aux   target system temperature for auxiliary dipoles
// uaux         auxiliary induced dipole value at each site
// upaux        auxiliary shadow induced dipoles at each site
// vaux         auxiliary induced dipole velocity at each site
// vpaux        auxiliary shadow dipole velocity at each site
// aaux         auxiliary induced dipole acceleration at each site
// apaux        auxiliary shadow dipole acceleration at each site
// use_ielscf   flag to use inertial extended Lagrangian method


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN int nfree_aux;
QCMD_EXTERN double tautemp_aux;
QCMD_EXTERN double kelvin_aux;
QCMD_EXTERN std::vector<std::vector<double>> uaux;
QCMD_EXTERN std::vector<std::vector<double>> upaux;
QCMD_EXTERN std::vector<std::vector<double>> vaux;
QCMD_EXTERN std::vector<std::vector<double>> vpaux;
QCMD_EXTERN std::vector<std::vector<double>> aaux;
QCMD_EXTERN std::vector<std::vector<double>> apaux;
QCMD_EXTERN bool use_ielscf;
