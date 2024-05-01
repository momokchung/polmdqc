// Author: Moses KJ Chung
// Year:   2024

#include "kdtree.h"
#include <algorithm>

namespace polmdqc
{
////////////////////////////////////////////////////////////
//                                                        //
//  kdtree  --  sort atoms geometrically using a kd-tree  //
//                                                        //
////////////////////////////////////////////////////////////

// "kdtree" consists of a set of procedures for sorting point
// geometrically using a kd-tree for AlphaMol2

struct alfatoms_cmp
{
    alfatoms_cmp (int index) : index_(index) {}

    bool operator() (const AlfAtom& atm1, const AlfAtom& atm2) const
    {
        return atm1.coord[index_] < atm2.coord[index_];
    }

    int index_;
};

void splitKDTree(std::vector<AlfAtom>& alfatoms, int begin, int end, int index, int nsplit, int nsplit_tot, std::vector<int>& v)
{
    if (nsplit >= nsplit_tot) return;

    const int dimensions = 3;
    int n = begin + (end - begin)/2;
    v.push_back(n);
    auto i = alfatoms.begin();
    std::nth_element(i + begin, i + n, i + end, alfatoms_cmp(index));

    index = (index+1) % dimensions;
    nsplit++;
    splitKDTree(alfatoms, begin, n, index, nsplit, nsplit_tot, v);
    splitKDTree(alfatoms, n , end, index, nsplit, nsplit_tot, v);
}

void kdTree(std::vector<AlfAtom>& alfatoms, int nsplit_tot, std::vector<int>& Nval)
{
    int begin = 0;
    int end = alfatoms.size();

    int nsplit= 0;
    int index = 0;

    std::vector<int> v;
    v.push_back(begin);
    v.push_back(end);

    splitKDTree(alfatoms, begin, end, index, nsplit, nsplit_tot, v);

    std::sort(v.begin(), v.end());

    for (int i = 0; i < v.size(); i++) Nval[i] = v[i];
}
}
