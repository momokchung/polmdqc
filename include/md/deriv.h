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

MDQC_EXTERN MDQCArray2D<real,3> desum;
MDQC_EXTERN MDQCArray2D<real,3> deb;
MDQC_EXTERN MDQCArray2D<real,3> dea;
MDQC_EXTERN MDQCArray2D<real,3> deba;
MDQC_EXTERN MDQCArray2D<real,3> deub;
MDQC_EXTERN MDQCArray2D<real,3> deaa;
MDQC_EXTERN MDQCArray2D<real,3> deopb;
MDQC_EXTERN MDQCArray2D<real,3> deopd;
MDQC_EXTERN MDQCArray2D<real,3> deid;
MDQC_EXTERN MDQCArray2D<real,3> deit;
MDQC_EXTERN MDQCArray2D<real,3> det;
MDQC_EXTERN MDQCArray2D<real,3> dept;
MDQC_EXTERN MDQCArray2D<real,3> debt;
MDQC_EXTERN MDQCArray2D<real,3> deat;
MDQC_EXTERN MDQCArray2D<real,3> dett;
MDQC_EXTERN MDQCArray2D<real,3> dev;
MDQC_EXTERN MDQCArray2D<real,3> der;
MDQC_EXTERN MDQCArray2D<real,3> dedsp;
MDQC_EXTERN MDQCArray2D<real,3> dec;
MDQC_EXTERN MDQCArray2D<real,3> decd;
MDQC_EXTERN MDQCArray2D<real,3> ded;
MDQC_EXTERN MDQCArray2D<real,3> dem;
MDQC_EXTERN MDQCArray2D<real,3> dep;
MDQC_EXTERN MDQCArray2D<real,3> dect;
MDQC_EXTERN MDQCArray2D<real,3> derxf;
MDQC_EXTERN MDQCArray2D<real,3> des;
MDQC_EXTERN MDQCArray2D<real,3> delf;
MDQC_EXTERN MDQCArray2D<real,3> deg;
MDQC_EXTERN MDQCArray2D<real,3> dex;

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

MDQC_EXTERN MDQCArray2D<real,3> ndesum;
MDQC_EXTERN MDQCArray2D<real,3> ndeb;
MDQC_EXTERN MDQCArray2D<real,3> ndea;
MDQC_EXTERN MDQCArray2D<real,3> ndeba;
MDQC_EXTERN MDQCArray2D<real,3> ndeub;
MDQC_EXTERN MDQCArray2D<real,3> ndeaa;
MDQC_EXTERN MDQCArray2D<real,3> ndeopb;
MDQC_EXTERN MDQCArray2D<real,3> ndeopd;
MDQC_EXTERN MDQCArray2D<real,3> ndeid;
MDQC_EXTERN MDQCArray2D<real,3> ndeit;
MDQC_EXTERN MDQCArray2D<real,3> ndet;
MDQC_EXTERN MDQCArray2D<real,3> ndept;
MDQC_EXTERN MDQCArray2D<real,3> ndebt;
MDQC_EXTERN MDQCArray2D<real,3> ndeat;
MDQC_EXTERN MDQCArray2D<real,3> ndett;
MDQC_EXTERN MDQCArray2D<real,3> ndev;
MDQC_EXTERN MDQCArray2D<real,3> nder;
MDQC_EXTERN MDQCArray2D<real,3> ndedsp;
MDQC_EXTERN MDQCArray2D<real,3> ndec;
MDQC_EXTERN MDQCArray2D<real,3> ndecd;
MDQC_EXTERN MDQCArray2D<real,3> nded;
MDQC_EXTERN MDQCArray2D<real,3> ndem;
MDQC_EXTERN MDQCArray2D<real,3> ndep;
MDQC_EXTERN MDQCArray2D<real,3> ndect;
MDQC_EXTERN MDQCArray2D<real,3> nderxf;
MDQC_EXTERN MDQCArray2D<real,3> ndes;
MDQC_EXTERN MDQCArray2D<real,3> ndelf;
MDQC_EXTERN MDQCArray2D<real,3> ndeg;
MDQC_EXTERN MDQCArray2D<real,3> ndex;
}
