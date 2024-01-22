// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "sizes.h"

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
MDQC_EXTERN MDQCArray<real> rsolv;
MDQC_EXTERN MDQCArray<real> rdescr;
MDQC_EXTERN MDQCArray<real> asolv;
MDQC_EXTERN MDQCArray<real> rborn;
MDQC_EXTERN MDQCArray<real> drb;
MDQC_EXTERN MDQCArray<real> drbp;
MDQC_EXTERN MDQCArray<real> drobc;
MDQC_EXTERN MDQCArray<real> gpol;
MDQC_EXTERN MDQCArray<real> shct;
MDQC_EXTERN MDQCArray<real> aobc;
MDQC_EXTERN MDQCArray<real> bobc;
MDQC_EXTERN MDQCArray<real> gobc;
MDQC_EXTERN MDQCArray<real> vsolv;
MDQC_EXTERN MDQCArray2D<real,maxclass> wace;
MDQC_EXTERN MDQCArray2D<real,maxclass> s2ace;
MDQC_EXTERN MDQCArray2D<real,maxclass> uace;
}
