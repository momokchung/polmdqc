// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  deriv  --  Cartesian coord derivative components  //
//                                                    //
////////////////////////////////////////////////////////

// desum   total energy Cartesian coordinate derivatives
// deb     bond stretch Cartesian coordinate derivatives
// dea     angle bend Cartesian coordinate derivatives
// deba    stretch-bend Cartesian coordinate derivatives
// deub    Urey-Bradley Cartesian coordinate derivatives
// deaa    angle-angle Cartesian coordinate derivatives
// deopb   out-of-plane bend Cartesian coordinate derivatives
// deopd   out-of-plane distance Cartesian coordinate derivatives
// deid    improper dihedral Cartesian coordinate derivatives
// deit    improper torsion Cartesian coordinate derivatives
// det     torsional Cartesian coordinate derivatives
// dept    pi-system torsion Cartesian coordinate derivatives
// debt    stretch-torsion Cartesian coordinate derivatives
// deat    angle-torsion Cartesian coordinate derivatives
// dett    torsion-torsion Cartesian coordinate derivatives
// dev     van der Waals Cartesian coordinate derivatives
// der     Pauli repulsion Cartesian coordinate derivatives
// dedsp   damped dispersion Cartesian coordinate derivatives
// dec     charge-charge Cartesian coordinate derivatives
// decd    charge-dipole Cartesian coordinate derivatives
// ded     dipole-dipole Cartesian coordinate derivatives
// dem     multipole Cartesian coordinate derivatives
// dep     polarization Cartesian coordinate derivatives
// dect    charge transfer Cartesian coordinate derivatives
// derxf   reaction field Cartesian coordinate derivatives
// des     solvation Cartesian coordinate derivatives
// delf    metal ligand field Cartesian coordinate derivatives
// deg     geometric restraint Cartesian coordinate derivatives
// dex     extra energy term Cartesian coordinate derivatives

MDQC_EXTERN std::vector<double> desum;
MDQC_EXTERN std::vector<double> deb;
MDQC_EXTERN std::vector<double> dea;
MDQC_EXTERN std::vector<double> deba;
MDQC_EXTERN std::vector<double> deub;
MDQC_EXTERN std::vector<double> deaa;
MDQC_EXTERN std::vector<double> deopb;
MDQC_EXTERN std::vector<double> deopd;
MDQC_EXTERN std::vector<double> deid;
MDQC_EXTERN std::vector<double> deit;
MDQC_EXTERN std::vector<double> det;
MDQC_EXTERN std::vector<double> dept;
MDQC_EXTERN std::vector<double> debt;
MDQC_EXTERN std::vector<double> deat;
MDQC_EXTERN std::vector<double> dett;
MDQC_EXTERN std::vector<double> dev;
MDQC_EXTERN std::vector<double> der;
MDQC_EXTERN std::vector<double> dedsp;
MDQC_EXTERN std::vector<double> dec;
MDQC_EXTERN std::vector<double> decd;
MDQC_EXTERN std::vector<double> ded;
MDQC_EXTERN std::vector<double> dem;
MDQC_EXTERN std::vector<double> dep;
MDQC_EXTERN std::vector<double> dect;
MDQC_EXTERN std::vector<double> derxf;
MDQC_EXTERN std::vector<double> des;
MDQC_EXTERN std::vector<double> delf;
MDQC_EXTERN std::vector<double> deg;
MDQC_EXTERN std::vector<double> dex;
}
