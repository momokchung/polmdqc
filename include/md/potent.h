// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  potent  --  usage of potential energy components  //
//                                                    //
////////////////////////////////////////////////////////

// use_bond    logical flag governing use of bond stretch potential
// use_angle   logical flag governing use of angle bend potential
// use_strbnd  logical flag governing use of stretch-bend potential
// use_urey    logical flag governing use of Urey-Bradley potential
// use_angang  logical flag governing use of angle-angle cross term
// use_opbend  logical flag governing use of out-of-plane bend term
// use_opdist  logical flag governing use of out-of-plane distance
// use_improp  logical flag governing use of improper dihedral term
// use_imptor  logical flag governing use of improper torsion term
// use_tors    logical flag governing use of torsional potential
// use_pitors  logical flag governing use of pi-system torsion term
// use_strtor  logical flag governing use of stretch-torsion term
// use_angtor  logical flag governing use of angle-torsion term
// use_tortor  logical flag governing use of torsion-torsion term
// use_vdw     logical flag governing use of van der Waals potential
// use_repel   logical flag governing use of Pauli repulsion term
// use_disp    logical flag governing use of dispersion potential
// use_charge  logical flag governing use of charge-charge potential
// use_chgdpl  logical flag governing use of charge-dipole potential
// use_dipole  logical flag governing use of dipole-dipole potential
// use_mpole   logical flag governing use of multipole potential
// use_polar   logical flag governing use of polarization term
// use_chgtrn  logical flag governing use of charge transfer term
// use_chgflx  logical flag governing use of charge flux term
// use_rxnfld  logical flag governing use of reaction field term
// use_solv    logical flag governing use of continuum solvation term
// use_metal   logical flag governing use of ligand field term
// use_geom    logical flag governing use of geometric restraints
// use_extra   logical flag governing use of extra potential term
// use_born    logical flag governing use of Born radii values
// use_orbit   logical flag governing use of pisystem computation
// use_mutate  logical flag governing use of hybrid potential terms

MDQC_EXTERN bool use_bond,use_angle;
MDQC_EXTERN bool use_strbnd,use_urey;
MDQC_EXTERN bool use_angang,use_opbend;
MDQC_EXTERN bool use_opdist,use_improp;
MDQC_EXTERN bool use_imptor,use_tors;
MDQC_EXTERN bool use_pitors,use_strtor;
MDQC_EXTERN bool use_angtor,use_tortor;
MDQC_EXTERN bool use_vdw,use_repel;
MDQC_EXTERN bool use_disp,use_charge;
MDQC_EXTERN bool use_chgdpl,use_dipole;
MDQC_EXTERN bool use_mpole,use_polar;
MDQC_EXTERN bool use_chgtrn,use_chgflx;
MDQC_EXTERN bool use_rxnfld,use_solv;
MDQC_EXTERN bool use_metal,use_geom;
MDQC_EXTERN bool use_extra,use_born;
MDQC_EXTERN bool use_orbit,use_mutate;
}
