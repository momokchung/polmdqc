///////////////////////////////////////////////////////
//                                                   //
//  torpot.h  --  torsional functional form details  //
//                                                   //
///////////////////////////////////////////////////////

// idihunit  convert improper dihedral energy to kcal/mole
// itorunit  convert improper torsion amplitudes to kcal/mole
// torsunit  convert torsional parameter amplitudes to kcal/mole
// ptorunit  convert pi-system torsion energy to kcal/mole
// storunit  convert stretch-torsion energy to kcal/mole
// atorunit  convert angle-torsion energy to kcal/mole
// ttorunit  convert torsion-torsion energy to kcal/mole


#pragma once
#include "macro.h"

QCMD_EXTERN double idihunit;
QCMD_EXTERN double itorunit;
QCMD_EXTERN double torsunit;
QCMD_EXTERN double ptorunit;
QCMD_EXTERN double storunit;
QCMD_EXTERN double atorunit;
QCMD_EXTERN double ttorunit;
