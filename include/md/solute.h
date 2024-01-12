// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  solute  --  continuum solvation model parameters  //
//                                                    //
////////////////////////////////////////////////////////

// doffset   dielectric offset to continuum solvation atomic radii
// onipr     probe radius to use with onion Born radius method
// p1        single-atom scale factor for analytical Still radii
// p2        1-2 interaction scale factor for analytical Still radii
// p3        1-3 interaction scale factor for analytical Still radii
// p4        nonbonded scale factor for analytical Still radii
// p5        soft cutoff parameter for analytical Still radii
// rsolv     atomic radius of each atom for continuum solvation
// rdescr    atomic radius of each atom for descreening
// asolv     atomic surface area solvation parameters
// rborn     Born radius of each atom for GB/SA solvation
// drb       solvation derivatives with respect to Born radii
// drbp      GK polarization derivatives with respect to Born radii
// drobc     chain rule term for Onufriev-Bashford-Case radii
// gpol      polarization self-energy values for each atom
// shct      overlap scale factors for Hawkins-Cramer-Truhlar radii
// aobc      alpha values for Onufriev-Bashford-Case radii
// bobc      beta values for Onufriev-Bashford-Case radii
// gobc      gamma values for Onufriev-Bashford-Case radii
// vsolv     atomic volume of each atom for use with ACE
// wace      "omega" values for atom class pairs for use with ACE
// s2ace     "sigma^2" values for atom class pairs for use with ACE
// uace      "mu" values for atom class pairs for use with ACE

MDQC_EXTERN real doffset,onipr;
MDQC_EXTERN real p1,p2,p3,p4,p5;
MDQC_EXTERN std::vector<real> rsolv;
MDQC_EXTERN std::vector<real> rdescr;
MDQC_EXTERN std::vector<real> asolv;
MDQC_EXTERN std::vector<real> rborn;
MDQC_EXTERN std::vector<real> drb;
MDQC_EXTERN std::vector<real> drbp;
MDQC_EXTERN std::vector<real> drobc;
MDQC_EXTERN std::vector<real> gpol;
MDQC_EXTERN std::vector<real> shct;
MDQC_EXTERN std::vector<real> aobc;
MDQC_EXTERN std::vector<real> bobc;
MDQC_EXTERN std::vector<real> gobc;
MDQC_EXTERN std::vector<real> vsolv;
MDQC_EXTERN std::vector<std::vector<real>> wace;
MDQC_EXTERN std::vector<std::vector<real>> s2ace;
MDQC_EXTERN std::vector<std::vector<real>> uace;
}
