// Author: Moses KJ Chung
// Year:   2024

#include "atomid.h"
#include "atoms.h"
#include "clusterqm.h"
#include "groupqm.h"

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  clusterqm  --  set user-defined groups of atoms  //
//                                                   //
///////////////////////////////////////////////////////

// "clusterqm" gets the partitioning of the system into groups
// and stores a list of the group to which each atom belongs

void clusterqm()
{
    // perform dynamic allocation of some global arrays
    igrpq.allocate(n);
    grpqmass.allocate(n);

    // find the first and last atom in each group
    int k = 0;
    igrpq[0][0] = 0;
    for (int i = 1; i < n; i++) {
        int j = grpqlist[i];
        if (j != k) {
            igrpq[k][1] = i - 1;
            k = j;
            igrpq[k][0] = i;
        }
    }
    igrpq[ngrpq-1][1] = n-1;

    // compute the mass of each group and the total mass
    tgrpqmass = 0.;
    for (int i = 0; i < ngrpq; i++) {
        grpqmass[i] = 0.;
        for (int k = igrpq[i][0]; k <= igrpq[i][1]; k++) {
            grpqmass[i] += mass[k];
        }
        tgrpqmass += grpqmass[i];
    }
}
}
