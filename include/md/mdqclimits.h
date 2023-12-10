// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////////////
//                                                        //
//  mdqclimits  --  interaction taper & cutoff distances  //
//                                                        //
////////////////////////////////////////////////////////////

// vdwcut      cutoff distance for van der Waals interactions
// repcut      cutoff distance for Pauli repulsion interactions
// dispcut     cutoff distance for dispersion interactions
// chgcut      cutoff distance for charge-charge interactions
// dplcut      cutoff distance for dipole-dipole interactions
// mpolecut    cutoff distance for atomic multipole interactions
// ctrncut     cutoff distance for charge transfer interactions
// vdwtaper    distance at which van der Waals switching begins
// reptaper    distance at which Pauli repulsion switching begins
// disptaper   distance at which dispersion switching begins
// chgtaper    distance at which charge-charge switching begins
// dpltaper    distance at which dipole-dipole switching begins
// mpoletaper  distance at which atomic multipole switching begins
// ctrntaper   distance at which charge transfer switching begins
// ewaldcut    cutoff distance for real space Ewald electrostatics
// dewaldcut   cutoff distance for real space Ewald dispersion
// usolvcut    cutoff distance for dipole solver preconditioner
// use_ewald   logical flag governing use of electrostatic Ewald
// use_dewald  logical flag governing use of dispersion Ewald
// use_lights  logical flag governing use of method of lights
// use_list    logical flag governing use of any neighbor lists
// use_vlist   logical flag governing use of van der Waals list
// use_dlist   logical flag governing use of dispersion list
// use_clist   logical flag governing use of charge list
// use_mlist   logical flag governing use of multipole list
// use_ulist   logical flag governing use of preconditioner list

MDQC_EXTERN double vdwcut,repcut;
MDQC_EXTERN double dispcut,chgcut;
MDQC_EXTERN double dplcut,mpolecut;
MDQC_EXTERN double ctrncut;
MDQC_EXTERN double vdwtaper,reptaper;
MDQC_EXTERN double disptaper,chgtaper;
MDQC_EXTERN double dpltaper,mpoletaper;
MDQC_EXTERN double ctrntaper;
MDQC_EXTERN double ewaldcut,dewaldcut;
MDQC_EXTERN double usolvcut;
MDQC_EXTERN bool use_ewald,use_dewald;
MDQC_EXTERN bool use_lights,use_list;
MDQC_EXTERN bool use_vlist,use_dlist;
MDQC_EXTERN bool use_clist,use_mlist;
MDQC_EXTERN bool use_ulist;
}
