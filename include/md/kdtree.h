// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "alfatom.h"
#include <vector>


namespace polmdqc
{
////////////////////////////////////////////////////////////
//                                                        //
//  kdtree  --  sort atoms geometrically using a kd-tree  //
//                                                        //
////////////////////////////////////////////////////////////

void kdTree(std::vector<AlfAtom>& alfatoms, int nsplit_tot, std::vector<int>& Nval);
}
