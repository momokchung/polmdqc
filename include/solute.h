//////////////////////////////////////////////////////////
//                                                      //
//  solute.h  --  continuum solvation model parameters  //
//                                                      //
//////////////////////////////////////////////////////////

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


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN double doffset,onipr;
QCMD_EXTERN double p1,p2,p3,p4,p5;
QCMD_EXTERN std::vector<double> rsolv;
QCMD_EXTERN std::vector<double> rdescr;
QCMD_EXTERN std::vector<double> asolv;
QCMD_EXTERN std::vector<double> rborn;
QCMD_EXTERN std::vector<double> drb;
QCMD_EXTERN std::vector<double> drbp;
QCMD_EXTERN std::vector<double> drobc;
QCMD_EXTERN std::vector<double> gpol;
QCMD_EXTERN std::vector<double> shct;
QCMD_EXTERN std::vector<double> aobc;
QCMD_EXTERN std::vector<double> bobc;
QCMD_EXTERN std::vector<double> gobc;
QCMD_EXTERN std::vector<double> vsolv;
QCMD_EXTERN std::vector<std::vector<double>> wace;
QCMD_EXTERN std::vector<std::vector<double>> s2ace;
QCMD_EXTERN std::vector<std::vector<double>> uace;
