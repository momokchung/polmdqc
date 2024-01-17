// Author: Moses KJ Chung
// Year:   2023

#include "angbnd.h"
#include "angles.h"
#include "atmlst.h"
#include "atoms.h"
#include "couple.h"
#include "fatal.h"

namespace polmdqc
{
////////////////////////////////////////////////
//                                            //
//  angles  --  locate and store bond angles  //
//                                            //
////////////////////////////////////////////////

// "angles" finds the total number of bond angles and stores
// the atom numbers of the atoms defining each angle; for
// each angle to a trivalent central atom, the third bonded
// atom is stored for use in out-of-plane bending

void angles()
{
    int ia,ib,ic;
    int maxang;

    // perform dynamic allocation of some global arrays
    maxang = 6 * n;
    if (iang.size() != 0) iang.resize(0);
    if (anglist.size() != 0) anglist.resize(0);
    if (balist.size() != 0) balist.resize(0);
    iang.resize(maxang, std::vector<int>(4));
    anglist.resize(n, std::vector<int>(maxval*(maxval-1)/2));
    balist.resize(maxang, std::vector<int>(2));

    // loop over all atoms, storing the atoms in each bond angle
    nangle = 0;
    for (int i = 0; i < n; i++) {
        int m = 0;
        for (int j = 0; j < n12[i]-1; j++) {
            for (int k = j+1; k < n12[i]; k++) {
                nangle++;
                if (nangle > maxang) {
                    printf("\n ANGLES  --  Too many Bond Angles; Increase MAXANG\n");
                    fatal();
                }
                anglist[i][m] = nangle;
                iang[nangle-1][0] = i12[i][j];
                iang[nangle-1][1] = i;
                iang[nangle-1][2] = i12[i][k];
                iang[nangle-1][3] = -1;
                m++;
            }
        }

        // set the out-of-plane atom for angles at trivalent centers
        if (n12[i] == 3) {
            iang[nangle-1][3] = i12[i][0];
            iang[nangle-2][3] = i12[i][1];
            iang[nangle-3][3] = i12[i][2];
        }
    }

    // store the numbers of the bonds comprising each bond angle
    for (int i = 0; i < nangle; i++) {
        int ia = iang[i][0];
        int ib = iang[i][1];
        int ic = iang[i][2];
        for (int k = 0; k < n12[ib]; k++) {
            if (i12[ib][k] == ia) balist[i][0] = bndlist[ib][k];
            if (i12[ib][k] == ic) balist[i][1] = bndlist[ib][k];
        }
    }
}
}
