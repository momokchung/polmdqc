// Author: Moses KJ Chung
// Year:   2024

#pragma once

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  delflip  --  flips used in Delaunay triangulation  //
//                                                     //
/////////////////////////////////////////////////////////

void flip_1_4(int ipoint, int itetra, int& tetra_last);
void flip();
}
