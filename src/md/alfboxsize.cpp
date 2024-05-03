// Author: Moses KJ Chung
// Year:   2024

#include "alfboxsize.h"

namespace polmdqc
{
/////////////////////////////////////////////////////////////
//                                                         //
//  alfboxsize  --  find smallest and largest coordinates  //
//                                                         //
/////////////////////////////////////////////////////////////

// "alfboxsize" finds the smallest and largest coordinates
// of the box

void alfboxsize(AlfAtom* alfatoms, int size, real& xmin, real& ymin, real& zmin, real& xmax, real& ymax, real& zmax, real& rmax)
{
    xmin = alfatoms[0].coord[0];
    xmax = alfatoms[0].coord[0];
    ymin = alfatoms[0].coord[1];
    ymax = alfatoms[0].coord[1];
    zmin = alfatoms[0].coord[2];
    zmax = alfatoms[0].coord[2];
    rmax = alfatoms[0].r;

    for (int i = 1; i < size; i++) {
        real x = alfatoms[i].coord[0];
        real y = alfatoms[i].coord[1];
        real z = alfatoms[i].coord[2];
        real r = alfatoms[i].r;

        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
        if (z < zmin) zmin = z;
        if (z > zmax) zmax = z;
        if (r > rmax) rmax = r;
    }
}
}
