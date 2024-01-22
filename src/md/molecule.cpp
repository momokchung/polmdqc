// Author: Moses KJ Chung
// Year:   2023

#include "atomid.h"
#include "atoms.h"
#include "couple.h"
#include "molcul.h"
#include "sort.h"
#include <algorithm>

namespace polmdqc
{
///////////////////////////////////////////////
//                                           //
//  molecule  --  assign atoms to molecules  //
//                                           //
///////////////////////////////////////////////

// "molecule" counts the molecules, assigns each atom to
// its molecule and computes the mass of each molecule

void molecule()
{
    std::vector<int> list;

    // perform dynamic allocation of some global arrays
    imol.allocate(n);
    kmol.allocate(n);
    molcule.allocate(n);
    molmass.allocate(n);

    // zero number of molecules and molecule membership list
    nmol = 0;
    for (int i = 0; i < n; i++) {
        molcule[i] = -1;
    }

    // assign each atom to its respective molecule
    for (int i = 0; i < n; i++) {
        if (molcule[i] == -1) {
            molcule[i] = nmol;
            nmol++;
        }
        int mi = molcule[i];
        for (int ii = 0; ii < n12[i]; ii++) {
            int j = i12[i][ii];
            int mj = molcule[j];
            if (mj == -1) {
                molcule[j] = mi;
            }
            else if (mi < mj) {
                nmol--;
                for (int k = 0; k < n; k++) {
                    int mk = molcule[k];
                    if (mk == mj) {
                        molcule[k] = mi;
                    }
                    else if (mk > mj) {
                        molcule[k] = mk - 1;
                    }
                }
            }
            else if (mi > mj) {
                nmol--;
                for (int k = 0; k < n; k++) {
                    int mk = molcule[k];
                    if (mk == mi) {
                        molcule[k] = mj;
                    }
                    else if (mk > mi) {
                        molcule[k] = mk - 1;
                    }
                }
                mi = mj;
            }
        }
    }

    // perform dynamic allocation of some local arrays
    list.resize(n);

    // pack atoms of each molecule into a contiguous indexed list
    for (int i = 0; i < n; i++) {
        list[i] = molcule[i];
    }
    sortKey(n, list, kmol.ptr());

    // find the first and last atom in each molecule
    int k = 0;
    imol[0][0] = 0;
    for (int i = 1; i < n; i++) {
        int j = list[i];
        if (j != k) {
            imol[k][1] = i - 1;
            k = j;
            imol[k][0] = i;
        }
    }
    imol[nmol-1][1] = n-1;

    // sort the list of atoms in each molecule by atom number
    for (int i = 0; i < nmol; i++) {
        k = imol[i][1] - imol[i][0] + 1;
        std::sort(kmol.ptr() + imol[i][0], kmol.ptr() + imol[i][0] + k);
    }

    // if all atomic masses are zero, set them all to unity
    for (int i = 0; i < n; i++) {
        if (mass[i] != 0.) goto label_10;
    }
    for (int i = 0; i < n; i++) {
        mass[i] = 1.;
    }
    label_10:

    // compute the mass of each molecule and the total mass
    totmass = 0.;
    for (int i = 0; i < nmol; i++) {
        molmass[i] = 0.;
        for (int k = imol[i][0]; k <= imol[i][1]; k++) {
            molmass[i] += mass[kmol[k]];
        }
        totmass += molmass[i];
    }
}
}
