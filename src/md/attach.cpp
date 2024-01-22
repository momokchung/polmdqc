// Author: Moses KJ Chung
// Year:   2023

#include "atoms.h"
#include "attach.h"
#include "couple.h"
#include "fatal.h"
#include "sort.h"

namespace polmdqc
{
////////////////////////////////////////////////
//                                            //
//  attach  --  setup of connectivity arrays  //
//                                            //
////////////////////////////////////////////////

// "attach" generates lists of 1-3, 1-4 and 1-5 connectivities
// starting from the previously determined list of attached
// atoms (ie, 1-2 connectivity)

void attach()
{
    // perform dynamic allocation of some global arrays
    int maxn13 = 3 * maxval;
    int maxn14 = 9 * maxval;
    int maxn15 = 27 * maxval;
    n13.allocate(n);
    n14.allocate(n);
    n15.allocate(n);
    i13.allocate(n);
    i14.allocate(n);
    i15.allocate(n);

    // loop over all atoms finding all the 1-3 relationships;
    // note "n12" and "i12" have already been setup elsewhere
    for (int i = 0; i < n; i++) {
        n13[i] = 0;
        std::vector<int> i13tmp;
        for (int j = 0; j < n12[i]; j++) {
            int jj = i12[i][j];
            for (int k = 0; k < n12[jj]; k++) {
                int kk = i12[jj][k];
                if (kk == i) goto label_10;
                for (int m = 0; m < n12[i]; m++){
                    if (kk == i12[i][m]) goto label_10;
                }
                i13tmp.push_back(kk);
                label_10:
                continue;
            }
        }
        int n13tmp = i13tmp.size();
        if (n13tmp > maxn13) {
            printf("\n ATTACH  --  Too many 1-3 Connected Atoms Attached to Atom%6d\n", i);
            fatal();
        }
        sortUnique(n13tmp, i13tmp, 0);
        n13[i] = n13tmp;
        for (int j = 0; j < n13tmp; j++){
            i13[i][j] = i13tmp[j];
        }
    }

    // loop over all atoms finding all the 1-4 relationships
    for (int i = 0; i < n; i++) {
        n14[i] = 0;
        std::vector<int> i14tmp;
        for (int j = 0; j < n13[i]; j++) {
            int jj = i13[i][j];
            for (int k = 0; k < n12[jj]; k++) {
                int kk = i12[jj][k];
                if (kk == i) goto label_30;
                for (int m = 0; m < n12[i]; m++){
                    if (kk == i12[i][m]) goto label_30;
                }
                for (int m = 0; m < n13[i]; m++){
                    if (kk == i13[i][m]) goto label_30;
                }
                i14tmp.push_back(kk);
                label_30:
                continue;
            }
        }
        int n14tmp = i14tmp.size();
        if (n14tmp > maxn14) {
            printf("\n ATTACH  --  Too many 1-4 Connected Atoms Attached to Atom%6d\n", i);
            fatal();
        }
        sortUnique(n14tmp, i14tmp, 0);
        n14[i] = n14tmp;
        for (int j = 0; j < n14tmp; j++){
            i14[i][j] = i14tmp[j];
        }
    }

    // loop over all atoms finding all the 1-5 relationships
    for (int i = 0; i < n; i++) {
        n15[i] = 0;
        std::vector<int> i15tmp;
        for (int j = 0; j < n14[i]; j++) {
            int jj = i14[i][j];
            for (int k = 0; k < n12[jj]; k++) {
                int kk = i12[jj][k];
                if (kk == i) goto label_50;
                for (int m = 0; m < n12[i]; m++){
                    if (kk == i12[i][m]) goto label_50;
                }
                for (int m = 0; m < n13[i]; m++){
                    if (kk == i13[i][m]) goto label_50;
                }
                for (int m = 0; m < n14[i]; m++){
                    if (kk == i14[i][m]) goto label_50;
                }
                i15tmp.push_back(kk);
                label_50:
                continue;
            }
        }
        int n15tmp = i15tmp.size();
        if (n15tmp > maxn15) {
            printf("\n ATTACH  --  Too many 1-5 Connected Atoms Attached to Atom%6d\n", i);
            fatal();
        }
        sortUnique(n15tmp, i15tmp, 0);
        n15[i] = n15tmp;
        for (int j = 0; j < n15tmp; j++){
            i15[i][j] = i15tmp[j];
        }
    }
}
}
