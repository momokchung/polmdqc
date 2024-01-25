// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "precision.h"

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
// te      torque on each atom (resolves gradient due to torque)

MDQC_EXTERN MDQCArray<real> desum;
MDQC_EXTERN MDQCArray<real> deb;
MDQC_EXTERN MDQCArray<real> dea;
MDQC_EXTERN MDQCArray<real> deba;
MDQC_EXTERN MDQCArray<real> deub;
MDQC_EXTERN MDQCArray<real> deaa;
MDQC_EXTERN MDQCArray<real> deopb;
MDQC_EXTERN MDQCArray<real> deopd;
MDQC_EXTERN MDQCArray<real> deid;
MDQC_EXTERN MDQCArray<real> deit;
MDQC_EXTERN MDQCArray<real> det;
MDQC_EXTERN MDQCArray<real> dept;
MDQC_EXTERN MDQCArray<real> debt;
MDQC_EXTERN MDQCArray<real> deat;
MDQC_EXTERN MDQCArray<real> dett;
MDQC_EXTERN MDQCArray<real> dev;
MDQC_EXTERN MDQCArray<real> der;
MDQC_EXTERN MDQCArray<real> dedsp;
MDQC_EXTERN MDQCArray<real> dec;
MDQC_EXTERN MDQCArray<real> decd;
MDQC_EXTERN MDQCArray<real> ded;
MDQC_EXTERN MDQCArray<real> dem;
MDQC_EXTERN MDQCArray<real> dep;
MDQC_EXTERN MDQCArray<real> dect;
MDQC_EXTERN MDQCArray<real> derxf;
MDQC_EXTERN MDQCArray<real> des;
MDQC_EXTERN MDQCArray<real> delf;
MDQC_EXTERN MDQCArray<real> deg;
MDQC_EXTERN MDQCArray<real> dex;
MDQC_EXTERN MDQCArray<real> te;

// ndesum   total energy Cartesian coordinate numerical derivatives
// ndeb     bond stretch Cartesian coordinate numerical derivatives
// ndea     angle bend Cartesian coordinate numerical derivatives
// ndeba    stretch-bend Cartesian coordinate numerical derivatives
// ndeub    Urey-Bradley Cartesian coordinate numerical derivatives
// ndeaa    angle-angle Cartesian coordinate numerical derivatives
// ndeopb   out-of-plane bend Cartesian coordinate numerical derivatives
// ndeopd   out-of-plane distance Cartesian coordinate numerical derivatives
// ndeid    improper dihedral Cartesian coordinate numerical derivatives
// ndeit    improper torsion Cartesian coordinate numerical derivatives
// ndet     torsional Cartesian coordinate numerical derivatives
// ndept    pi-system torsion Cartesian coordinate numerical derivatives
// ndebt    stretch-torsion Cartesian coordinate numerical derivatives
// ndeat    angle-torsion Cartesian coordinate numerical derivatives
// ndett    torsion-torsion Cartesian coordinate numerical derivatives
// ndev     van der Waals Cartesian coordinate numerical derivatives
// nder     Pauli repulsion Cartesian coordinate numerical derivatives
// ndedsp   damped dispersion Cartesian coordinate numerical derivatives
// ndec     charge-charge Cartesian coordinate numerical derivatives
// ndecd    charge-dipole Cartesian coordinate numerical derivatives
// nded     dipole-dipole Cartesian coordinate numerical derivatives
// ndem     multipole Cartesian coordinate numerical derivatives
// ndep     polarization Cartesian coordinate numerical derivatives
// ndect    charge transfer Cartesian coordinate numerical derivatives
// nderxf   reaction field Cartesian coordinate numerical derivatives
// ndes     solvation Cartesian coordinate numerical derivatives
// ndelf    metal ligand field Cartesian coordinate numerical derivatives
// ndeg     geometric restraint Cartesian coordinate numerical derivatives
// ndex     extra energy term Cartesian coordinate numerical derivatives

MDQC_EXTERN MDQCArray<real> ndesum;
MDQC_EXTERN MDQCArray<real> ndeb;
MDQC_EXTERN MDQCArray<real> ndea;
MDQC_EXTERN MDQCArray<real> ndeba;
MDQC_EXTERN MDQCArray<real> ndeub;
MDQC_EXTERN MDQCArray<real> ndeaa;
MDQC_EXTERN MDQCArray<real> ndeopb;
MDQC_EXTERN MDQCArray<real> ndeopd;
MDQC_EXTERN MDQCArray<real> ndeid;
MDQC_EXTERN MDQCArray<real> ndeit;
MDQC_EXTERN MDQCArray<real> ndet;
MDQC_EXTERN MDQCArray<real> ndept;
MDQC_EXTERN MDQCArray<real> ndebt;
MDQC_EXTERN MDQCArray<real> ndeat;
MDQC_EXTERN MDQCArray<real> ndett;
MDQC_EXTERN MDQCArray<real> ndev;
MDQC_EXTERN MDQCArray<real> nder;
MDQC_EXTERN MDQCArray<real> ndedsp;
MDQC_EXTERN MDQCArray<real> ndec;
MDQC_EXTERN MDQCArray<real> ndecd;
MDQC_EXTERN MDQCArray<real> nded;
MDQC_EXTERN MDQCArray<real> ndem;
MDQC_EXTERN MDQCArray<real> ndep;
MDQC_EXTERN MDQCArray<real> ndect;
MDQC_EXTERN MDQCArray<real> nderxf;
MDQC_EXTERN MDQCArray<real> ndes;
MDQC_EXTERN MDQCArray<real> ndelf;
MDQC_EXTERN MDQCArray<real> ndeg;
MDQC_EXTERN MDQCArray<real> ndex;
}
