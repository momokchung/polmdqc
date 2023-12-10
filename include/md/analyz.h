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

MDQC_EXTERN std::vector<double> aesum;
MDQC_EXTERN std::vector<double> aeb;
MDQC_EXTERN std::vector<double> aea;
MDQC_EXTERN std::vector<double> aeba;
MDQC_EXTERN std::vector<double> aeub;
MDQC_EXTERN std::vector<double> aeaa;
MDQC_EXTERN std::vector<double> aeopb;
MDQC_EXTERN std::vector<double> aeopd;
MDQC_EXTERN std::vector<double> aeid;
MDQC_EXTERN std::vector<double> aeit;
MDQC_EXTERN std::vector<double> aet;
MDQC_EXTERN std::vector<double> aept;
MDQC_EXTERN std::vector<double> aebt;
MDQC_EXTERN std::vector<double> aeat;
MDQC_EXTERN std::vector<double> aett;
MDQC_EXTERN std::vector<double> aev;
MDQC_EXTERN std::vector<double> aer;
MDQC_EXTERN std::vector<double> aedsp;
MDQC_EXTERN std::vector<double> aec;
MDQC_EXTERN std::vector<double> aecd;
MDQC_EXTERN std::vector<double> aed;
MDQC_EXTERN std::vector<double> aem;
MDQC_EXTERN std::vector<double> aep;
MDQC_EXTERN std::vector<double> aect;
MDQC_EXTERN std::vector<double> aerxf;
MDQC_EXTERN std::vector<double> aes;
MDQC_EXTERN std::vector<double> aelf;
MDQC_EXTERN std::vector<double> aeg;
MDQC_EXTERN std::vector<double> aex;
}
