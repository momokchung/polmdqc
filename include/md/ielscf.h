// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  ielscf  --  extended Lagrangian induced dipoles  //
//                                                   //
///////////////////////////////////////////////////////

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

MDQC_EXTERN int nfree_aux;
MDQC_EXTERN real tautemp_aux;
MDQC_EXTERN real kelvin_aux;
MDQC_EXTERN std::vector<std::vector<real>> uaux;
MDQC_EXTERN std::vector<std::vector<real>> upaux;
MDQC_EXTERN std::vector<std::vector<real>> vaux;
MDQC_EXTERN std::vector<std::vector<real>> vpaux;
MDQC_EXTERN std::vector<std::vector<real>> aaux;
MDQC_EXTERN std::vector<std::vector<real>> apaux;
MDQC_EXTERN bool use_ielscf;
}
