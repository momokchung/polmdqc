// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  torpot  --  torsional functional form details  //
//                                                 //
/////////////////////////////////////////////////////

// idihunit  convert improper dihedral energy to kcal/mole
// itorunit  convert improper torsion amplitudes to kcal/mole
// torsunit  convert torsional parameter amplitudes to kcal/mole
// ptorunit  convert pi-system torsion energy to kcal/mole
// storunit  convert stretch-torsion energy to kcal/mole
// atorunit  convert angle-torsion energy to kcal/mole
// ttorunit  convert torsion-torsion energy to kcal/mole

MDQC_EXTERN real idihunit;
MDQC_EXTERN real itorunit;
MDQC_EXTERN real torsunit;
MDQC_EXTERN real ptorunit;
MDQC_EXTERN real storunit;
MDQC_EXTERN real atorunit;
MDQC_EXTERN real ttorunit;
}
