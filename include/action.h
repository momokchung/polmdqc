///////////////////////////////////////////////////////////
//                                                       //
//  action.h  --  total number of each energy term type  //
//                                                       //
///////////////////////////////////////////////////////////

// neb     number of bond stretch energy terms computed
// nea     number of angle bend energy terms computed
// neba    number of stretch-bend energy terms computed
// neub    number of Urey-Bradley energy terms computed
// neaa    number of angle-angle energy terms computed
// neopb   number of out-of-plane bend energy terms computed
// neopd   number of out-of-plane distance energy terms computed
// neid    number of improper dihedral energy terms computed
// neit    number of improper torsion energy terms computed
// net     number of torsional energy terms computed
// nept    number of pi-system torsion energy terms computed
// nebt    number of stretch-torsion energy terms computed
// neat    number of angle-torsion energy terms computed
// nett    number of torsion-torsion energy terms computed
// nev     number of van der Waals energy terms computed
// ner     number of Pauli repulsion energy terms computed
// nedsp   number of dispersion energy terms computed
// nec     number of charge-charge energy terms computed
// necd    number of charge-dipole energy terms computed
// ned     number of dipole-dipole energy terms computed
// nem     number of multipole energy terms computed
// nep     number of polarization energy terms computed
// nect    number of charge transfer energy terms computed
// news    number of Ewald summation energy terms computed
// nerxf   number of reaction field energy terms computed
// nes     number of solvation energy terms computed
// nelf    number of metal ligand field energy terms computed
// neg     number of geometric restraint energy terms computed
// nex     number of extra energy terms computed


#pragma once
#include "macro.h"

QCMD_EXTERN int neb,nea,neba,neub;
QCMD_EXTERN int neaa,neopb,neopd;
QCMD_EXTERN int neid,neit,net,nept;
QCMD_EXTERN int nebt,neat,nett,nev;
QCMD_EXTERN int ner,nedsp,nec,necd;
QCMD_EXTERN int ned,nem,nep,nect;
QCMD_EXTERN int news,nerxf,nes,nelf;
QCMD_EXTERN int neg,nex;
