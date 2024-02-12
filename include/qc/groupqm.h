// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////////////
//                                                        //
//  groupqm  --  partitioning of system into atom groups  //
//                                                        //
////////////////////////////////////////////////////////////

// ngrpq       total number of atom groups in the system
// tgrpqmass   total mass of all the atoms
// grpqlist    number of the group to which each atom belongs
// igrpq       first and last atom of each group in the list
// grpqchg     molecular charge for each molecule in the system
// grpqmult    molecular multiplicity for each molecule in the system
// grpqmass    total mass of all the atoms in each group

MDQC_EXTERN int ngrpq;
MDQC_EXTERN realq tgrpqmass;
MDQC_EXTERN MDQCArray<int> grpqlist;
MDQC_EXTERN MDQCArray2D<int,2> igrpq;
MDQC_EXTERN MDQCArray<int> grpqchg;
MDQC_EXTERN MDQCArray<int> grpqmult;
MDQC_EXTERN MDQCArray<real> grpqmass;
}
