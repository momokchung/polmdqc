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

MDQC_EXTERN double idihunit;
MDQC_EXTERN double itorunit;
MDQC_EXTERN double torsunit;
MDQC_EXTERN double ptorunit;
MDQC_EXTERN double storunit;
MDQC_EXTERN double atorunit;
MDQC_EXTERN double ttorunit;
}
