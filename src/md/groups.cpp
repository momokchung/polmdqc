// Author: Moses KJ Chung
// Year:   2023

#include "group.h"
#include "groups.h"
#include <algorithm>
#include <cmath>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  groups  --  group membership of set of atoms  //
//                                                //
////////////////////////////////////////////////////

// "groups" tests a set of atoms to see if all are members of a
// single atom group or a pair of atom groups; if so, then the
// correct intra- or intergroup weight is assigned
// 
// note the default group-based interaction weight is 1.0; only
// interactions involving two or fewer groups can be scaled

void groups(bool& proceed, real& weigh, int ia, int ib, int ic, int id, int ie, int ig)
{
    int iga,igb,igc;
    int igd,ige,igg;
    int nset;
    int gmax,gmin;

    // determine the number of atoms in the set to be compared
    nset = 0;
    weigh = 1.;
    if (ig != -1) nset = 6;
    else if (ie != -1) nset = 5;
    else if (id != -1) nset = 4;
    else if (ic != -1) nset = 3;
    else if (ib != -1) nset = 2;
    else if (ia != -1) nset = 1;

    // check group membership for a set containing one atom
    if (nset == 1) {
        iga = grplist[ia]+1;
        weigh = wgrp[iga][iga];
    }

    // check group membership for a set containing two atoms
    else if (nset == 2) {
        iga = grplist[ia]+1;
        igb = grplist[ib]+1;
        weigh = wgrp[igb][iga];
    }

    // check group membership for a set containing three atoms
    else if (nset == 3) {
        iga = grplist[ia]+1;
        igb = grplist[ib]+1;
        igc = grplist[ic]+1;
        if (iga==igb or igb==igc) {
            weigh = wgrp[igc][iga];
        }
        else if (iga == igc) {
            weigh = wgrp[igb][iga];
        }
    }

    // check group membership for a set containing four atoms
    else if (nset == 4) {
        iga = grplist[ia]+1;
        igb = grplist[ib]+1;
        igc = grplist[ic]+1;
        igd = grplist[id]+1;
        gmin = std::min({iga,igb,igc,igd});
        gmax = std::max({iga,igb,igc,igd});
        if ((iga==gmin or iga==gmax) and
            (igb==gmin or igb==gmax) and
            (igc==gmin or igc==gmax) and
            (igd==gmin or igd==gmax)) weigh = wgrp[gmax][gmin];
    }

    // check group membership for a set containing five atoms
    else if (nset == 5) {
        iga = grplist[ia]+1;
        igb = grplist[ib]+1;
        igc = grplist[ic]+1;
        igd = grplist[id]+1;
        ige = grplist[ie]+1;
        gmin = std::min({iga,igb,igc,igd,ige});
        gmax = std::max({iga,igb,igc,igd,ige});
        if ((iga==gmin or iga==gmax) and
            (igb==gmin or igb==gmax) and
            (igc==gmin or igc==gmax) and
            (igd==gmin or igd==gmax) and
            (ige==gmin or ige==gmax)) weigh = wgrp[gmax][gmin];
    }

    // check group membership for a set containing five atoms
    else if (nset == 6) {
        iga = grplist[ia]+1;
        igb = grplist[ib]+1;
        igc = grplist[ic]+1;
        igd = grplist[id]+1;
        ige = grplist[ie]+1;
        igg = grplist[ig]+1;
        gmin = std::min({iga,igb,igc,igd,ige,igg});
        gmax = std::max({iga,igb,igc,igd,ige,igg});
        if ((iga==gmin or iga==gmax) and
            (igb==gmin or igb==gmax) and
            (igc==gmin or igc==gmax) and
            (igd==gmin or igd==gmax) and
            (ige==gmin or ige==gmax) and
            (igg==gmin or igg==gmax))  weigh = wgrp[gmax][gmin];
    }

    // interaction will be used if its group has nonzero weight
    if (weigh == 0.) proceed = false;
    else proceed = true;
}
}
