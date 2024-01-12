// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

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

MDQC_EXTERN std::vector<real> aesum;
MDQC_EXTERN std::vector<real> aeb;
MDQC_EXTERN std::vector<real> aea;
MDQC_EXTERN std::vector<real> aeba;
MDQC_EXTERN std::vector<real> aeub;
MDQC_EXTERN std::vector<real> aeaa;
MDQC_EXTERN std::vector<real> aeopb;
MDQC_EXTERN std::vector<real> aeopd;
MDQC_EXTERN std::vector<real> aeid;
MDQC_EXTERN std::vector<real> aeit;
MDQC_EXTERN std::vector<real> aet;
MDQC_EXTERN std::vector<real> aept;
MDQC_EXTERN std::vector<real> aebt;
MDQC_EXTERN std::vector<real> aeat;
MDQC_EXTERN std::vector<real> aett;
MDQC_EXTERN std::vector<real> aev;
MDQC_EXTERN std::vector<real> aer;
MDQC_EXTERN std::vector<real> aedsp;
MDQC_EXTERN std::vector<real> aec;
MDQC_EXTERN std::vector<real> aecd;
MDQC_EXTERN std::vector<real> aed;
MDQC_EXTERN std::vector<real> aem;
MDQC_EXTERN std::vector<real> aep;
MDQC_EXTERN std::vector<real> aect;
MDQC_EXTERN std::vector<real> aerxf;
MDQC_EXTERN std::vector<real> aes;
MDQC_EXTERN std::vector<real> aelf;
MDQC_EXTERN std::vector<real> aeg;
MDQC_EXTERN std::vector<real> aex;
}
