// Author: Moses KJ Chung
// Year:   2024

#include "alforder.h"
#include "alfp.h"
#include "dlauny2.h"
#include "grid.h"
#include "hilbert.h"
#include "kdtree.h"
#include <cmath>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  alforder  --  reorder/partition atoms for AlphaMol2  //
//                                                       //
///////////////////////////////////////////////////////////

// "alforder" reorders/partitions atoms based on Sort3D,
// BRIO, Split, or KDTree algorithms

void alforder(real xmin, real ymin, real zmin, real xmax, real ymax, real zmax, real rmax, int nthreads, std::vector<int>& Nval)
{
    if (alfsort == AlfSort::None) return;

    int depth = 0;

    if (alfsort == AlfSort::Sort3D) {
        sort3DHilbert(&alfatoms[0], alfatoms.size(), 0, 0, xmin, xmax, ymin, ymax, zmin, zmax, depth);
    }
    else if (alfsort == AlfSort::BRIO) {
        brioHilbert(&alfatoms[0], alfatoms.size(), xmin, xmax, ymin, ymax, zmin, zmax, depth);
    }
    else if (alfsort == AlfSort::Split) {
        splitGrid(&alfatoms[0], alfatoms.size(), xmin, xmax, ymin, ymax, zmin, zmax, nthreads, Nval);
    }
    else if (alfsort == AlfSort::KDTree) {
        int nsplit = (int) std::log2((double) nthreads);
        if (nsplit > 0) kdTree(alfatoms, nsplit, Nval);
    }
}
}
