// Author: Moses KJ Chung
// Year:   2024

#include "angbnd.h"
#include "angles.h"
#include "atoms.h"
#include "bitor.h"
#include "bitors.h"
#include "bndstr.h"
#include "bonds.h"
#include "couple.h"
#include "inform.h"
#include "fatal.h"
#include "ring.h"
#include "rings.h"
#include "tors.h"
#include "torsions.h"
#include <algorithm>
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////
//                                           //
//  rings  --  locate and store small rings  //
//                                           //
///////////////////////////////////////////////

// "rings" searches the structure for small rings and stores
// their constituent atoms, and optionally reduces large rings
// into their component smaller rings

// note by default reducible rings are not removed since they
// are needed for force field parameter assignment

void rings()
{
    int imax;
    int ia,ib,ic,id,ie,ig,ih;
    int list1,list2,list3,list4;
    int maxring;
    std::vector<int> list;
    bool reduce;

    // zero out the number of small rings in the structure
    reduce = false;
    nring3 = 0;
    nring4 = 0;
    nring5 = 0;
    nring6 = 0;
    nring7 = 0;

    // parse to find bonds, angles, torsions and bitorsions
    if (nbond == 0) bonds();
    if (nangle == 0) angles();
    if (ntors == 0) torsions();
    if (nbitor == 0) bitors();

    // perform dynamic allocation of some global arrays
    maxring = 10000;
    iring3.allocate(maxring);
    iring4.allocate(maxring);
    iring5.allocate(maxring);
    iring6.allocate(maxring);
    iring7.allocate(maxring);

    // search for and store all of the 3-membered rings
    for (int i = 0; i < nangle; i++) {
        ia = iang[i][0];
        ib = iang[i][1];
        ic = iang[i][2];
        if (ib<ia and ib<ic) {
            for (int j = 0; j < n12[ia]; j++) {
                if (i12[ia][j] == ic) {
                    if (nring3 >= maxring) {
                        printf("\n RINGS  --  Too many 3-Membered Rings; Increase MAXRING\n");
                        fatal();
                    }
                    iring3[nring3][0] = ia;
                    iring3[nring3][1] = ib;
                    iring3[nring3][2] = ic;
                    nring3 = nring3 + 1;
                    continue;
                }
            }
        }
    }

    // perform dynamic allocation of some local arrays
    list.resize(n);

    // search for and store all of the 4-membered rings
    for (int i = 0; i < n; i++) {
        list[i] = 0;
    }
    for (int i = 0; i < ntors; i++) {
        ia = itors[i][0];
        ib = itors[i][1];
        ic = itors[i][2];
        id = itors[i][3];
        if (ia<ic and id<ib) {
            for (int j = 0; j < n12[ia]; j++) {
                if (i12[ia][j] == id) {
                    if (nring4 >= maxring) {
                        printf("\n RINGS  --  Too many 4-Membered Rings; Increase MAXRING\n");
                        fatal();
                    }
                    iring4[nring4][0] = ia;
                    iring4[nring4][1] = ib;
                    iring4[nring4][2] = ic;
                    iring4[nring4][3] = id;
                    nring4 += 1;

                    // remove the ring if it is reducible into smaller rings
                    if (reduce) {
                        list[ia] = nring4;
                        list[ib] = nring4;
                        list[ic] = nring4;
                        list[id] = nring4;
                        for (int m = 0; m < nring3; m++) {
                            list1 = list[iring3[m][0]];
                            list2 = list[iring3[m][1]];
                            list3 = list[iring3[m][2]];
                            if (list1==nring4 and list2==nring4 and list3==nring4) {
                                nring4 -= 1;
                                list[ia] = 0;
                                list[ib] = 0;
                                list[ic] = 0;
                                list[id] = 0;
                                goto label_40;
                            }
                        }
                    }
                    goto label_40;
                }
            }
            label_40:;
        }
    }

    // search for and store all of the 5-membered rings
    for (int i = 0; i < n; i++) {
        list[i] = 0;
    }
    for (int i = 0; i < nbitor; i++) {
        ia = ibitor[i][0];
        ib = ibitor[i][1];
        ic = ibitor[i][2];
        id = ibitor[i][3];
        ie = ibitor[i][4];
        if (ia<id and ie<ib and std::min(ia,ie)<ic) {
            for (int j = 0; j < n12[ia]; j++) {
                if (i12[ia][j] == ie) {
                    if (nring5 >= maxring) {
                        printf("\n RINGS  --  Too many 5-Membered Rings; Increase MAXRING\n");
                        fatal();
                    }
                    iring5[nring5][0] = ia;
                    iring5[nring5][1] = ib;
                    iring5[nring5][2] = ic;
                    iring5[nring5][3] = id;
                    iring5[nring5][4] = ie;
                    nring5 += 1;

                    // remove the ring if it is reducible into smaller rings
                    if (reduce) {
                        list[ia] = nring5;
                        list[ib] = nring5;
                        list[ic] = nring5;
                        list[id] = nring5;
                        list[ie] = nring5;
                        for (int m = 0; m < nring3; m++) {
                            list1 = list[iring3[m][0]];
                            list2 = list[iring3[m][1]];
                            list3 = list[iring3[m][2]];
                            if (list1==nring5 and list2==nring5 and list3==nring5) {
                                nring5 -= 1;
                                list[ia] = 0;
                                list[ib] = 0;
                                list[ic] = 0;
                                list[id] = 0;
                                list[ie] = 0;
                                goto label_60;
                            }
                        }
                    }
                    goto label_60;
                }
            }
            label_60:;
        }
    }

    // search for and store all of the 6-membered rings
    for (int i = 0; i < n; i++) {
        list[i] = 0;
    }
    for (int i = 0; i < nbitor; i++) {
        ia = ibitor[i][0];
        ib = ibitor[i][1];
        ic = ibitor[i][2];
        id = ibitor[i][3];
        ie = ibitor[i][4];
        imax = std::max({ia,ib,ic,id,ie});
        for (int j = 0; j < n12[ia]; j++) {
            ig = i12[ia][j];
            if (ig > imax) {
                for (int k = 0; k < n12[ie]; k++) {
                    if (i12[ie][k] == ig) {
                        if (nring6 >= maxring) {
                            printf("\n RINGS  --  Too many 6-Membered Rings; Increase MAXRING\n");
                            fatal();
                        }
                        iring6[nring6][0] = ia;
                        iring6[nring6][1] = ib;
                        iring6[nring6][2] = ic;
                        iring6[nring6][3] = id;
                        iring6[nring6][4] = ie;
                        iring6[nring6][5] = ig;
                        nring6 += 1;

                        // remove the ring if it is reducible into smaller rings
                        if (reduce) {
                            list[ia] = nring6;
                            list[ib] = nring6;
                            list[ic] = nring6;
                            list[id] = nring6;
                            list[ie] = nring6;
                            list[ig] = nring6;
                            for (int m = 0; m < nring3; m++) {
                                list1 = list[iring3[m][0]];
                                list2 = list[iring3[m][1]];
                                list3 = list[iring3[m][2]];
                                if (list1==nring6 and list2==nring6 and list3==nring6) {
                                    nring6 -= 1;
                                    list[ia] = 0;
                                    list[ib] = 0;
                                    list[ic] = 0;
                                    list[id] = 0;
                                    list[ie] = 0;
                                    list[ig] = 0;
                                    goto label_80;
                                }
                            }
                            for (int m = 0; m < nring4; m++) {
                                list1 = list[iring4[m][0]];
                                list2 = list[iring4[m][1]];
                                list3 = list[iring4[m][2]];
                                list4 = list[iring4[m][3]];
                                if (list1==nring6 and list2==nring6 and list3==nring6 and list4==nring6) {
                                    nring6 -= 1;
                                    list[ia] = 0;
                                    list[ib] = 0;
                                    list[ic] = 0;
                                    list[id] = 0;
                                    list[ie] = 0;
                                    list[ig] = 0;
                                    goto label_80;
                                }
                            }
                        }
                        label_80:;
                    }
                }
            }
        }
    }

    // search for and store all of the 7-membered rings
    for (int i = 0; i < n; i++) {
        list[i] = 0;
    }
    for (int i = 0; i < nbitor; i++) {
        ia = ibitor[i][0];
        ib = ibitor[i][1];
        ic = ibitor[i][2];
        id = ibitor[i][3];
        ie = ibitor[i][4];
        imax = std::max({ia,ib,ic,id,ie});
        for (int j = 0; j < n12[ia]; j++) {
            ih = i12[ia][j];
            for (int k = 0; k < n12[ie]; k++) {
                ig = i12[ie][k];
                if (ig!=id and ih!=ib and ((ig>imax and ih>ie) or (ih>imax and ig>ia))) {
                    for (int kk = 0; kk < n12[ig]; kk++) {
                        if (i12[ig][kk] == ih) {
                            if (nring7 >= maxring) {
                                printf("\n RINGS  --  Too many 7-Membered Rings; Increase MAXRING\n");
                                fatal();
                            }
                            iring7[nring7][0] = ia;
                            iring7[nring7][1] = ib;
                            iring7[nring7][2] = ic;
                            iring7[nring7][3] = id;
                            iring7[nring7][4] = ie;
                            iring7[nring7][5] = ig;
                            iring7[nring7][6] = ih;
                            nring7 += 1;

                            // remove the ring if it is reducible into smaller rings
                            if (reduce) {
                                list[ia] = nring7;
                                list[ib] = nring7;
                                list[ic] = nring7;
                                list[id] = nring7;
                                list[ie] = nring7;
                                list[ig] = nring7;
                                list[ih] = nring7;
                                for (int m = 0; m < nring3; m++) {
                                    list1 = list[iring3[m][0]];
                                    list2 = list[iring3[m][1]];
                                    list3 = list[iring3[m][2]];
                                    if (list1==nring7 and list2==nring7 and list3==nring7) {
                                        nring7 -= 1;
                                        list[ia] = 0;
                                        list[ib] = 0;
                                        list[ic] = 0;
                                        list[id] = 0;
                                        list[ie] = 0;
                                        list[ig] = 0;
                                        list[ih] = 0;
                                        goto label_100;
                                    }
                                }
                                for (int m = 0; m < nring4; m++) {
                                    list1 = list[iring4[m][0]];
                                    list2 = list[iring4[m][1]];
                                    list3 = list[iring4[m][2]];
                                    list4 = list[iring4[m][3]];
                                    if (list1==nring7 and list2==nring7 and list3==nring7 and list4==nring7) {
                                        nring7 -=  1;
                                        list[ia] = 0;
                                        list[ib] = 0;
                                        list[ic] = 0;
                                        list[id] = 0;
                                        list[ie] = 0;
                                        list[ig] = 0;
                                        list[ih] = 0;
                                        goto label_100;
                                    }
                                }
                            }
                            label_100:;
                        }
                    }
                }
            }
        }
    }

    // print out lists of the small rings in the structure
    if (debug) {
        int space1,space2;
        if (nring3 > 0) {
            space1 = 3;
            space2 = 14;
            printf("\n Three-Membered Rings in the Structure :");
            printf("\n\n%*sRing%*sAtoms in Ring\n\n",space1,"",space2,"");
            for (int i = 0; i < nring3; i++) {
                printf("%6d       %7d%7d%7d\n",i+1,iring3[i][0]+1,iring3[i][1]+1,iring3[i][2]+1);
            }
        }
        if (nring4 > 0) {
            space1 = 3;
            space2 = 17;
            printf("\n Four-Membered Rings in the Structure :");
            printf("\n\n%*sRing%*sAtoms in Ring\n\n",space1,"",space2,"");
            for (int i = 0; i < nring4; i++) {
                printf("%6d       %7d%7d%7d%7d\n",i+1,iring4[i][0]+1,iring4[i][1]+1,iring4[i][2]+1,iring4[i][3]+1);
            }
        }
        if (nring5 > 0) {
            space1 = 3;
            space2 = 20;
            printf("\n Five-Membered Rings in the Structure :");
            printf("\n\n%*sRing%*sAtoms in Ring\n\n",space1,"",space2,"");
            for (int i = 0; i < nring5; i++) {
                printf("%6d       %7d%7d%7d%7d%7d\n",i+1,iring5[i][0]+1,iring5[i][1]+1,iring5[i][2]+1,iring5[i][3]+1,iring5[i][4]+1);
            }
        }
        if (nring6 > 0) {
            space1 = 3;
            space2 = 23;
            printf("\n Six-Membered Rings in the Structure :");
            printf("\n\n%*sRing%*sAtoms in Ring\n\n",space1,"",space2,"");
            for (int i = 0; i < nring6; i++) {
                printf("%6d       %7d%7d%7d%7d%7d%7d\n",i+1,iring6[i][0]+1,iring6[i][1]+1,iring6[i][2]+1,iring6[i][3]+1,iring6[i][4]+1,iring6[i][5]+1);
            }
        }
        if (nring7 > 0) {
            space1 = 3;
            space2 = 26;
            printf("\n Seven-Membered Rings in the Structure :");
            printf("\n\n%*sRing%*sAtoms in Ring\n\n",space1,"",space2,"");
            for (int i = 0; i < nring7; i++) {
                printf("%6d       %7d%7d%7d%7d%7d%7d%7d\n",i+1,iring7[i][0]+1,iring7[i][1]+1,iring7[i][2]+1,iring7[i][3]+1,iring7[i][4]+1,iring7[i][5]+1,iring7[i][6]+1);
            }
        }
    }
}
}