// Author: Moses KJ Chung
// Year:   2023

#include "atomid.h"
#include "atoms.h"
#include "bound.h"
#include "cluster.h"
#include "cutoffs.h"
#include "fatal.h"
#include "getnumb.h"
#include "gettext.h"
#include "group.h"
#include "inform.h"
#include "keys.h"
#include "molcul.h"
#include "mdqclimits.h"
#include "sort.h"
#include "upcase.h"
#include <algorithm>
#include <cmath>
#include <sstream>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  cluster  --  set user-defined groups of atoms  //
//                                                 //
/////////////////////////////////////////////////////

// "cluster" gets the partitioning of the system into groups
// and stores a list of the group to which each atom belongs

void cluster()
{
    int next,size;
    int gnum,ga,gb;
    std::vector<int> list;
    real wg;
    bool header;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // perform dynamic allocation of some global arrays
    if (igrp.size() == 0) igrp.resize(maxgrp+1, std::vector<int>(2)); // same
    if (grpmass.size() == 0) grpmass.resize(maxgrp+1);
    if (wgrp.size() == 0) wgrp.resize(maxgrp+1, std::vector<real>(maxgrp+1, 1.)); // same
    if (kgrp.size() != 0) (kgrp.resize(0)); // same
    if (grplist.size() != 0) (grplist.resize(0)); // same
    kgrp.resize(n, -1);
    grplist.resize(n, -1);

    // set defaults for the group atom list and weight options
    use_group = false;
    use_intra = false;
    use_inter = false;
    ngrp = 0;
    for (int i = 0; i <= maxgrp; i++) {
        igrp[i][0] = 0;
        igrp[i][1] = -1;
    }

    // perform dynamic allocation of some local arrays
    size = std::max(100,n);
    list.resize(size);

    // get any keywords containing atom group definitions
    for (int j = 0; j < nkey; j++) {
        next = 0;
        record = keyline[j];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "GROUP") {
            use_group = true;
            gnum = 0;
            for (int i = 0; i < size; i++) {
               list[i] = 0;
            }
            getnumb(record,gnum,next);
            if (gnum > maxgrp) {
                printf("\n CLUSTER  --  Too many Atom Groups; Increase MAXGRP\n");
                fatal();
            }
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            int listInt;
            int counter = 0;
            while (iss>>listInt and counter<size) {
                list[counter] = listInt;
                counter++;
            }
            int i = 0;
            while (list[i] != 0) {
                if (list[i] > 0) {
                    grplist[list[i]-1] = gnum-1;
                    i += 1;
                }
                else {
                    for (int k = std::abs(list[i]); k <= std::abs(list[i+1]); k++) {
                        grplist[k-1] = gnum-1;
                    }
                    i += 2;
                }
            }
        }

        // get any keywords with weights for group interactions
        else if (keyword == "GROUP-MOLECULE") {
            use_group = true;
            use_inter = true;
            use_intra = false;
            if (nmol > maxgrp) {
                printf("\n CLUSTER  --  Too many Atom Groups; Increase MAXGRP\n");
                fatal();
            }
            for (int i = 0; i < nmol; i++) {
                for (int k = imol[i][0]; k <= imol[i][1]; k++) {
                    grplist[kmol[k]] = i;
                }
            }
        }

        // get any keywords with weights for group interactions
        else if (keyword == "GROUP-SELECT") {
            ga = 0;
            gb = 0;
            wg = -1.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ga >> gb >> wg;
            if (wg < 0.)  wg = 1.;
            wgrp[gb][ga] = wg;
            wgrp[ga][gb] = wg;
            use_inter = false;
        }

        // get keywords to select common sets of group interactions
        else if (keyword == "GROUP-INTRA") {
            use_intra = true;
            use_inter = false;
        }
        else if (keyword == "GROUP-INTER") {
            use_inter = true;
            use_intra = false;
        }
    }

    // pack atoms of each group into a contiguous indexed list
    if (use_group) {
        for (int i = 0; i < n; i++) {
            list[i] = grplist[i];
        }
        sortKey(n, list, kgrp);

        // find the first and last atom in each of the groups
        int k = list[0]+1;
        igrp[k][0] = 0;
        int j;
        for (int i = 0; i < n; i++) {
            j = list[i]+1;
            if (j != k) {
                igrp[k][1] = i - 1;
                igrp[j][0] = i;
                k = j;
            }
            ngrp = std::max(ngrp, j);
        }
        igrp[j][1] = n-1;

        // sort the list of atoms in each group by atom number
        for (int i = 0; i <= ngrp; i++) {
            int size = igrp[i][1] - igrp[i][0] + 1;
            if (igrp[i][0] != -1) std::sort(kgrp.begin() + igrp[i][0], kgrp.begin() + igrp[i][0] + size);
        }
    }

    // use only intragroup or intergroup interactions if selected
    if (use_intra) {
        for (int i = 0; i <= ngrp; i++) {
            for (int j = 0; j <= ngrp; j++) {
                wgrp[i][j] = 0.;
            }
            wgrp[i][i] = 1.;
        }
    }
    if (use_inter) {
        for (int i = 0; i <= ngrp; i++) {
            for (int j = 0; j <= ngrp; j++) {
                wgrp[i][j] = 1.;
            }
            wgrp[i][i] = 0.;
        }
    }

    // disable consideration of interactions with any empty groups
    for (int i = 0; i <= ngrp; i++) {
        size = igrp[i][1] - igrp[i][0] + 1;
        if (size == 0) {
            for (int j = 0; j <= ngrp; j++) {
                wgrp[i][j] = 0.;
                wgrp[j][i] = 0.;
            }
        }
    }

    // turn off bounds and replicas for intragroup calculations
    if (use_intra) {
        use_bounds = false;
        use_replica = false;
        cutoffs();
    }

    // compute the total mass of all atoms in each group
    for (int i = 1; i <= ngrp; i++) {
        grpmass[i] = 0.;
        for (int j = igrp[i][0]; j <= igrp[i][1]; j++) {
            grpmass[i] += mass[kgrp[j]];
        }
    }

    // output the final list of atoms in each group
    if (use_group and debug) {
        for (int i = 1; i <= ngrp; i++) {
            size = igrp[i][1] - igrp[i][0] + 1;
            if (size != 0) {
                printf("\n List of Atoms in Group%3d :\n\n   ", i);
                int counter = 0;
                for (int j = igrp[i][0]; j <= igrp[i][1]; j++) {
                    printf("%7d", kgrp[j]+1);
                    counter++;
                    if ((counter) % 10 == 0) printf("\n   ");
                }
                printf("\n");
            }
        }
    }

    // output the weights for intragroup and intergroup interactions
    if (use_group and debug) {
        header = true;
        for (int i = 0; i <= ngrp; i++) {
            for (int j = i; j <= ngrp; j++) {
                if (wgrp[i][j] != 0.) {
                    if (header) {
                        header = false;
                        printf("\n Active Sets of Intra- and InterGroup Interactions :\n\n");
                        printf("           Groups               Type              Weight\n\n");
                    }
                    if (i == j) {
                        printf("     %6d%6d            IntraGroup     %12.4f\n", i, j, wgrp[i][j]);
                    }
                    else {
                        printf("     %6d%6d            InterGroup     %12.4f\n", i, j, wgrp[i][j]);
                    }
                }
            }
        }
    }
}
}
