// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  energi  --  individual potential energy components  //
//                                                      //
//////////////////////////////////////////////////////////

// esum   total potential energy of the system
// eb     bond stretch potential energy of the system
// ea     angle bend potential energy of the system
// eba    stretch-bend potential energy of the system
// eub    Urey-Bradley potential energy of the system
// eaa    angle-angle potential energy of the system
// eopb   out-of-plane bend potential energy of the system
// eopd   out-of-plane distance potential energy of the system
// eid    improper dihedral potential energy of the system
// eit    improper torsion potential energy of the system
// et     torsional potential energy of the system
// ept    pi-system torsion potential energy of the system
// ebt    stretch-torsion potential energy of the system
// eat    angle-torsion potential energy of the system
// ett    torsion-torsion potential energy of the system
// ev     van der Waals potential energy of the system
// er     Pauli repulsion potential energy of the system
// edsp   dampled dispersion potential energy of the system
// ec     charge-charge potential energy of the system
// ecd    charge-dipole potential energy of the system
// ed     dipole-dipole potential energy of the system
// em     atomic multipole potential energy of the system
// ep     polarization potential energy of the system
// ect    charge transfer potential energy of the system
// erxf   reaction field potential energy of the system
// es     solvation potential energy of the system
// elf    metal ligand field potential energy of the system
// eg     geometric restraint potential energy of the system
// ex     extra term potential energy of the system

MDQC_EXTERN double esum,eb,ea;
MDQC_EXTERN double eba,eub,eaa;
MDQC_EXTERN double eopb,eopd,eid;
MDQC_EXTERN double eit,et,ept;
MDQC_EXTERN double ebt,eat,ett;
MDQC_EXTERN double ev,er,edsp;
MDQC_EXTERN double ec,ecd,ed;
MDQC_EXTERN double em,ep,ect;
MDQC_EXTERN double erxf,es,elf;
MDQC_EXTERN double eg,ex;
}
