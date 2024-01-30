// Author: Moses KJ Chung
// Year:   2024

#include "chkring.h"
#include "couple.h"

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  chkring  --  check atom set for small rings  //
//                                               //
///////////////////////////////////////////////////

// "chkring" tests an atom or a set of connected atoms for
// their presence within a single 3- to 6-membered ring

void chkring(int& iring, int ia, int ib, int ic, int id)
{
    int i,m,q,r;
    int nset;

    // initialize the ring size and number of atoms to test
    iring = 0;
    nset = 0;
    if (ia >= 0) nset = 1;
    if (ib >= 0) nset = 2;
    if (ic >= 0) nset = 3;
    if (id >= 0) nset = 4;

    // cannot be in a ring if the terminal atoms are univalent
    if (nset == 1) {
        if (n12[ia] <= 1) nset = 0;
    }
    else if (nset == 2) {
        if (std::min(n12[ia],n12[ib]) <= 1) nset = 0;
    }
    else if (nset == 3) {
        if (std::min(n12[ia],n12[ic]) <= 1) nset = 0;
    }
    else if (nset == 4) {
        if (std::min(n12[ia],n12[id]) <= 1) nset = 0;
    }

    // check the input atoms for sequential connectivity
    if (nset > 1) {
        for (int j = 0; j < n12[ia]; j++) {
            i = i12[ia][j];
            if (ib == i) {
                if (nset == 2) goto label_10;
                for (int k = 0; k < n12[ib]; k++) {
                    m = i12[ib][k];
                    if (ic == m) {
                        if (nset == 3) goto label_10;
                        for (int p = 0; p < n12[ic]; p++) {
                            q = i12[ic][p];
                            if (id == q) goto label_10;
                        }
                    }
                }
            }
        }
        nset = 0;
    }
    label_10:;

    // check for an atom contained inside a small ring
    if (nset == 1) {
        for (int j = 0; j < n12[ia]-1; j++) {
            i = i12[ia][j];
            for (int k = j+1; k < n12[ia]; k++) {
                m = i12[ia][k];
                for (int p = 0; p < n12[i]; p++) {
                    if (m == i12[i][p]) {
                        iring = 3;
                        goto label_end;
                    }
                }
            }
        }
        for (int j = 0; j < n12[ia]-1; j++) {
            i = i12[ia][j];
            for (int k = j+1; k < n12[ia]; k++) {
                m = i12[ia][k];
                for (int p = 0; p < n12[i]; p++) {
                    r = i12[i][p];
                    if (r != ia) {
                        for (int q = 0; q < n12[m]; q++) {
                            if (r == i12[m][q]) {
                                iring = 4;
                                goto label_end;
                            }
                        }
                    }
                }
            }
        }
        for (int j = 0; j < n13[ia]-1; j++) {
            i = i13[ia][j];
            for (int k = j+1; k < n13[ia]; k++) {
                m = i13[ia][k];
                for (int p = 0; p < n12[i]; p++) {
                    if (m == i12[i][p]) {
                        iring = 5;
                        goto label_end;
                    }
                }
                for (int p = 0; p < n13[i]; p++) {
                    if (m == i13[i][p]) {
                        iring = 6;
                        goto label_end;
                    }
                }
            }
        }
    }

    // check for a bond contained inside a small ring
    else if (nset == 2) {
        for (int j = 0; j < n12[ia]; j++) {
            i = i12[ia][j];
            for (int k = 0; k < n12[ib]; k++) {
                if (i == i12[ib][k]) {
                    iring = 3;
                    goto label_end;
                }
            }
        }
        for (int j = 0; j < n12[ia]; j++) {
            i = i12[ia][j];
            if (ib != i) {
                for (int k = 0; k < n12[ib]; k++) {
                    m = i12[ib][k];
                    if (ia != m) {
                        for (int p = 0; p < n12[i]; p++) {
                            if (m == i12[i][p]) {
                                iring = 4;
                                goto label_end;
                            }
                        }
                    }
                }
            }
        }
        for (int j = 0; j < n13[ia]; j++) {
            i = i13[ia][j];
            for (int k = 0; k < n13[ib]; k++) {
                if (i == i13[ib][k]) {
                    iring = 5;
                    goto label_end;
                }
            }
        }
        for (int j = 0; j < n12[ia]; j++) {
            i = i12[ia][j];
            if (ib != i) {
                for (int k = 0; k < n13[ib]; k++) {
                    m = i13[ib][k];
                    for (int p = 0; p < n13[i]; p++) {
                        if (m == i13[i][p]) {
                            iring = 6;
                            for (int q = 0; q < n12[ia]; q++) {
                                if (m == i12[ia][q]) iring = 0;
                            }
                            if (iring == 6) goto label_end;
                        }
                    }
                }
            }
        }
    }

    // check for an angle contained inside a small ring
    else if (nset == 3) {
        for (int j = 0; j < n12[ia]; j++) {
            if (ic == i12[ia][j]) {
                iring = 3;
                goto label_end;
            }
        }
        for (int j = 0; j < n12[ia]; j++) {
            i = i12[ia][j];
            if (ib != i) {
                for (int k = 0; k < n12[ic]; k++) {
                    if (i == i12[ic][k]) {
                        iring = 4;
                        goto label_end;
                    }
                }
            }
        }
        for (int j = 0; j < n12[ia]; j++) {
            i = i12[ia][j];
            if (ib != i) {
                for (int k = 0; k < n13[ic]; k++) {
                    if (i == i13[ic][k]) {
                        iring = 5;
                        goto label_end;
                    }
                }
            }
        }
        for (int j = 0; j < n13[ia]; j++) {
            i = i13[ia][j];
            if (ic != i) {
                for (int k = 0; k < n13[ic]; k++) {
                    if (i == i13[ic][k]) {
                        iring = 6;
                        goto label_end;
                    }
                }
            }
        }
    } 

    // check for a torsion contained inside a small ring
    else if (nset == 4) {
        for (int j = 0; j < n12[ia]; j++) {
            if (id == i12[ia][j]) {
                iring = 4;
                goto label_end;
            }
        }
        for (int j = 0; j < n12[ia]; j++) {
            i = i12[ia][j];
            if (ib != i) {
                for (int k = 0; k < n12[id]; k++) {
                    if (i == i12[id][k]) {
                        iring = 5;
                        goto label_end;
                    }
                }
            }
        }
        for (int j = 0; j < n12[ia]; j++) {
            i = i12[ia][j];
            if (ib != i) {
                for (int k = 0; k < n13[id]; k++) {
                    if (i == i13[id][k]) {
                        iring = 6;
                        goto label_end;
                    }
                }
            }
        }
    }
    label_end:;
}
}
