// Author: Moses KJ Chung
// Year:   2023

#include "atmlst.h"
#include "atoms.h"
#include "bndstr.h"
#include "bonds.h"
#include "couple.h"
#include "fatal.h"

namespace polmdqc
{
//////////////////////////////////////////////////
//                                              //
//  bonds  --  locate and store covalent bonds  //
//                                              //
//////////////////////////////////////////////////

// "bonds" finds the total number of covalent bonds and
// stores the atom numbers of the atoms defining each bond

void bonds()
{
    int maxbnd;

    // perform dynamic allocation of some global arrays
    maxbnd = 4 * n;
    ibnd.allocate(maxbnd);
    bndlist.allocate(n);

    // loop over all atoms, storing the atoms in each bond
    nbond = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n12[i]; j++) {
            int k = i12[i][j];
            if (i < k) {
                nbond++;
                if (nbond > maxbnd) {
                    printf("\n BONDS  --  Too many Bonds; Increase MAXBND\n");
                    fatal();
                }
                ibnd[nbond-1][0] = i;
                ibnd[nbond-1][1] = k;
                bndlist[i][j] = nbond;
                for (int m = 0; m < n12[k]; m++) {
                    if (i == i12[k][m]) {
                        bndlist[k][m] = nbond;
                        break;
                    }
                }
            }
        }
    }
}
}
