////////////////////////////////////////////////////////////
//                                                        //
//  analyz.h  --  energy components partitioned to atoms  //
//                                                        //
////////////////////////////////////////////////////////////

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


#pragma once
#include "macro.h"
#include <vector>

QCMD_EXTERN std::vector<double> aesum;
QCMD_EXTERN std::vector<double> aeb;
QCMD_EXTERN std::vector<double> aea;
QCMD_EXTERN std::vector<double> aeba;
QCMD_EXTERN std::vector<double> aeub;
QCMD_EXTERN std::vector<double> aeaa;
QCMD_EXTERN std::vector<double> aeopb;
QCMD_EXTERN std::vector<double> aeopd;
QCMD_EXTERN std::vector<double> aeid;
QCMD_EXTERN std::vector<double> aeit;
QCMD_EXTERN std::vector<double> aet;
QCMD_EXTERN std::vector<double> aept;
QCMD_EXTERN std::vector<double> aebt;
QCMD_EXTERN std::vector<double> aeat;
QCMD_EXTERN std::vector<double> aett;
QCMD_EXTERN std::vector<double> aev;
QCMD_EXTERN std::vector<double> aer;
QCMD_EXTERN std::vector<double> aedsp;
QCMD_EXTERN std::vector<double> aec;
QCMD_EXTERN std::vector<double> aecd;
QCMD_EXTERN std::vector<double> aed;
QCMD_EXTERN std::vector<double> aem;
QCMD_EXTERN std::vector<double> aep;
QCMD_EXTERN std::vector<double> aect;
QCMD_EXTERN std::vector<double> aerxf;
QCMD_EXTERN std::vector<double> aes;
QCMD_EXTERN std::vector<double> aelf;
QCMD_EXTERN std::vector<double> aeg;
QCMD_EXTERN std::vector<double> aex;
