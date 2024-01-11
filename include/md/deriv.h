// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "precision.h"
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

MDQC_EXTERN std::vector<std::vector<real>> desum;
MDQC_EXTERN std::vector<std::vector<real>> deb;
MDQC_EXTERN std::vector<std::vector<real>> dea;
MDQC_EXTERN std::vector<std::vector<real>> deba;
MDQC_EXTERN std::vector<std::vector<real>> deub;
MDQC_EXTERN std::vector<std::vector<real>> deaa;
MDQC_EXTERN std::vector<std::vector<real>> deopb;
MDQC_EXTERN std::vector<std::vector<real>> deopd;
MDQC_EXTERN std::vector<std::vector<real>> deid;
MDQC_EXTERN std::vector<std::vector<real>> deit;
MDQC_EXTERN std::vector<std::vector<real>> det;
MDQC_EXTERN std::vector<std::vector<real>> dept;
MDQC_EXTERN std::vector<std::vector<real>> debt;
MDQC_EXTERN std::vector<std::vector<real>> deat;
MDQC_EXTERN std::vector<std::vector<real>> dett;
MDQC_EXTERN std::vector<std::vector<real>> dev;
MDQC_EXTERN std::vector<std::vector<real>> der;
MDQC_EXTERN std::vector<std::vector<real>> dedsp;
MDQC_EXTERN std::vector<std::vector<real>> dec;
MDQC_EXTERN std::vector<std::vector<real>> decd;
MDQC_EXTERN std::vector<std::vector<real>> ded;
MDQC_EXTERN std::vector<std::vector<real>> dem;
MDQC_EXTERN std::vector<std::vector<real>> dep;
MDQC_EXTERN std::vector<std::vector<real>> dect;
MDQC_EXTERN std::vector<std::vector<real>> derxf;
MDQC_EXTERN std::vector<std::vector<real>> des;
MDQC_EXTERN std::vector<std::vector<real>> delf;
MDQC_EXTERN std::vector<std::vector<real>> deg;
MDQC_EXTERN std::vector<std::vector<real>> dex;
}
