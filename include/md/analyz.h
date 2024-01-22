// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  analyz  --  energy components partitioned to atoms  //
//                                                      //
//////////////////////////////////////////////////////////

// aesum   total potential energy partitioned over atoms
// aeb     bond stretch energy partitioned over atoms
// aea     angle bend energy partitioned over atoms
// aeba    stretch-bend energy partitioned over atoms
// aeub    Urey-Bradley energy partitioned over atoms
// aeaa    angle-angle energy partitioned over atoms
// aeopb   out-of-plane bend energy partitioned over atoms
// aeopd   out-of-plane distance energy partitioned over atoms
// aeid    improper dihedral energy partitioned over atoms
// aeit    improper torsion energy partitioned over atoms
// aet     torsional energy partitioned over atoms
// aept    pi-system torsion energy partitioned over atoms
// aebt    stretch-torsion energy partitioned over atoms
// aeat    angle-torsion energy partitioned over atoms
// aett    torsion-torsion energy partitioned over atoms
// aev     van der Waals energy partitioned over atoms
// aer     Pauli repulsion energy partitioned over atoms
// aedsp   damped dispersion energy partitioned over atoms
// aec     charge-charge energy partitioned over atoms
// aecd    charge-dipole energy partitioned over atoms
// aed     dipole-dipole energy partitioned over atoms
// aem     multipole energy partitioned over atoms
// aep     polarization energy partitioned over atoms
// aect    charge transfer energy partitioned over atoms
// aerxf   reaction field energy partitioned over atoms
// aes     solvation energy partitioned over atoms
// aelf    metal ligand field energy partitioned over atoms
// aeg     geometric restraint energy partitioned over atoms
// aex     extra energy term partitioned over atoms

MDQC_EXTERN MDQCArray<real> aesum;
MDQC_EXTERN MDQCArray<real> aeb;
MDQC_EXTERN MDQCArray<real> aea;
MDQC_EXTERN MDQCArray<real> aeba;
MDQC_EXTERN MDQCArray<real> aeub;
MDQC_EXTERN MDQCArray<real> aeaa;
MDQC_EXTERN MDQCArray<real> aeopb;
MDQC_EXTERN MDQCArray<real> aeopd;
MDQC_EXTERN MDQCArray<real> aeid;
MDQC_EXTERN MDQCArray<real> aeit;
MDQC_EXTERN MDQCArray<real> aet;
MDQC_EXTERN MDQCArray<real> aept;
MDQC_EXTERN MDQCArray<real> aebt;
MDQC_EXTERN MDQCArray<real> aeat;
MDQC_EXTERN MDQCArray<real> aett;
MDQC_EXTERN MDQCArray<real> aev;
MDQC_EXTERN MDQCArray<real> aer;
MDQC_EXTERN MDQCArray<real> aedsp;
MDQC_EXTERN MDQCArray<real> aec;
MDQC_EXTERN MDQCArray<real> aecd;
MDQC_EXTERN MDQCArray<real> aed;
MDQC_EXTERN MDQCArray<real> aem;
MDQC_EXTERN MDQCArray<real> aep;
MDQC_EXTERN MDQCArray<real> aect;
MDQC_EXTERN MDQCArray<real> aerxf;
MDQC_EXTERN MDQCArray<real> aes;
MDQC_EXTERN MDQCArray<real> aelf;
MDQC_EXTERN MDQCArray<real> aeg;
MDQC_EXTERN MDQCArray<real> aex;
}
