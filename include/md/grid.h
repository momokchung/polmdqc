// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  grid  --  sort atoms geometrically using 3D grid  //
//                                                    //
////////////////////////////////////////////////////////

void splitGrid(AlfAtom *alfatoms, int size, real xmin, real xmax, real ymin, real ymax, real zmin, real zmax, int ncube, std::vector<int>& Nval);
}
