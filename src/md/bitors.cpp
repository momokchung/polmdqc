// Author: Moses KJ Chung
// Year:   2023

#include "angbnd.h"
#include "atoms.h"
#include "bitor.h"
#include "bitors.h"
#include "couple.h"
#include "fatal.h"

namespace polmdqc
{
///////////////////////////////////////////////
//                                           //
//  bitors  --  locate and store bitorsions  //
//                                           //
///////////////////////////////////////////////

// "bitors" finds the total number of bitorsions as pairs
// of adjacent torsional angles, and the numbers of the five
// atoms defining each bitorsion

void bitors()
{
    int maxbitor;

    // perform dynamic allocation of some global arrays
    maxbitor = 54 * n;
    ibitor.allocate(maxbitor);

    // loop over all angles, storing the atoms in each bitorsion
    nbitor = 0;
    for (int i = 0; i < nangle; i++) {
        int ib = iang[i][0];
        int ic = iang[i][1];
        int id = iang[i][2];
        for (int j = 0; j < n12[ib]; j++) {
            int ia = i12[ib][j];
            if (ia!=ic and ia!=id) {
                for (int k = 0; k < n12[id]; k++) {
                    int ie = i12[id][k];
                    if (ie!=ic and ie!=ib and ie!=ia) {
                        nbitor++;
                        if (nbitor > maxbitor) {
                            printf("\n BITORS  --  Too many Adjacent Torsions; Increase MAXBITOR\n");
                            fatal();
                        }
                        ibitor[nbitor-1][0] = ia;
                        ibitor[nbitor-1][1] = ib;
                        ibitor[nbitor-1][2] = ic;
                        ibitor[nbitor-1][3] = id;
                        ibitor[nbitor-1][4] = ie;
                    }
                }
            }
        }
    }
}
}
