//////////////////////////////////////////////////////////////
//                                                          //
//  qcmdlimits.h  --  interaction taper & cutoff distances  //
//                                                          //
//////////////////////////////////////////////////////////////

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


#pragma once
#include "macro.h"

QCMD_EXTERN double vdwcut,repcut;
QCMD_EXTERN double dispcut,chgcut;
QCMD_EXTERN double dplcut,mpolecut;
QCMD_EXTERN double ctrncut;
QCMD_EXTERN double vdwtaper,reptaper;
QCMD_EXTERN double disptaper,chgtaper;
QCMD_EXTERN double dpltaper,mpoletaper;
QCMD_EXTERN double ctrntaper;
QCMD_EXTERN double ewaldcut,dewaldcut;
QCMD_EXTERN double usolvcut;
QCMD_EXTERN bool use_ewald,use_dewald;
QCMD_EXTERN bool use_lights,use_list;
QCMD_EXTERN bool use_vlist,use_dlist;
QCMD_EXTERN bool use_clist,use_mlist;
QCMD_EXTERN bool use_ulist;
