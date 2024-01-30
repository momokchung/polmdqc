// Author: Moses KJ Chung
// Year:   2023

#include "atoms.h"
#include "chgpen.h"
#include "chkpole.h"
#include "couple.h"
#include "expol.h"
#include "fatal.h"
#include "getnumb.h"
#include "gettext.h"
#include "inform.h"
#include "keys.h"
#include "kpolar.h"
#include "kpolpr.h"
#include "kpolr.h"
#include "mplpot.h"
#include "mpole.h"
#include "numeral.h"
#include "polar.h"
#include "polgrp.h"
#include "polopt.h"
#include "polpcg.h"
#include "polpot.h"
#include "poltcg.h"
#include "potent.h"
#include "sort.h"
#include "upcase.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  kpolar  --  assign polarizability parameters  //
//                                                //
////////////////////////////////////////////////////

// "kpolar" assigns atomic dipole polarizabilities to the atoms
// within the structure and processes any new or changed values
// 
// literature references:
// 
// A. C. Simmonett, F. C. Pickard IV, J. W. Ponder and B. R. Brooks,
// "An Empirical Extrapolation Scheme for Efficient Treatment of
// Induced Dipoles", Journal of Chemical Physics, 145, 164101 (2016)
// [OPT method]
// 
// F. Aviat, L. Lagardere and J.-P. Piquemal, "The Truncated
// Conjugate Gradient (TCG), a Non-Iterative/Fixed-Cost Strategy for
// Computing Polarization in Molecular Dynamics: Fast Evaluation of
// Analytical Forces", Journal of Chemical Physics, 147, 161724
// (2018) [TCG method]

void kpolar()
{
    int i,j,k;
    int ii,kk;
    int ia,ib,it;
    int next,size;
    int nlist,npg;
    int number;
    int pg[maxval];
    std::vector<int> list;
    std::vector<int> rlist;
    real pol,thl,thd;
    real sixth;
    bool header;
    std::string pa,pb;
    std::string blank,pt;
    std::string keyword;
    std::string text;
    std::string record;
    std::string string;
    std::istringstream iss;

    // set the default values for polarization variables
    polprt = false;

    // set defaults for PCG induced dipole parameters
    pcgprec = true;
    pcgguess = true;
    pcgpeek = 1.;

    // set defaults for TCG induced dipole parameters
    tcgorder = 0;
    tcgguess = true;
    tcgpeek = 1.;
    if (poltyp == "TCG") poltyp = "TCG2";
    if (poltyp == "TCG0") {
        poltyp = "DIRECT";
    }
    else if (poltyp == "TCG1") {
        poltyp = "TCG";
        tcgorder = 1;
    }
    else if (poltyp.substr(0,3) == "TCG") {
        poltyp = "TCG";
        tcgorder = 2;
    }

    // perform dynamic allocation of some global arrays
    copt.allocate(maxopt+1);
    copm.allocate(maxopt+1);

    // set defaults for OPT induced dipole coefficients
    optorder = 0;
    for (int i = 0; i < maxopt; i++) {
        copt[i] = 0.;
        copm[i] = 0.;
    }
    if (poltyp == "OPT") poltyp = "OPT4";
    if (poltyp == "OPT1") {
        copt[0] = 0.530;
        copt[1] = 0.604;
    }
    else if (poltyp == "OPT2") {
        copt[0] = 0.042;
        copt[1] = 0.635;
        copt[2] = 0.414;
    }
    else if (poltyp == "OPT3") {
        copt[0] = -0.132;
        copt[1] = 0.218;
        copt[2] = 0.637;
        copt[3] = 0.293;
    }
    else if (poltyp == "OPT4") {
        copt[0] = -0.071;
        copt[1] = -0.096;
        copt[2] = 0.358;
        copt[3] = 0.587;
        copt[4] = 0.216;
    }
    else if (poltyp == "OPT5") {
        copt[0] = -0.005;
        copt[1] = -0.129;
        copt[2] = -0.026;
        copt[3] = 0.465;
        copt[4] = 0.528;
        copt[5] = 0.161;
    }
    else if (poltyp == "OPT6") {
        copt[0] = 0.014;
        copt[1] = -0.041;
        copt[2] = -0.172;
        copt[3] = 0.073;
        copt[4] = 0.535;
        copt[5] = 0.467;
        copt[6] = 0.122;
    }

    // perform dynamic allocation of some local arrays
    list.resize(n);

    // set defaults for numbers and lists of polarizable atoms
    nlist = 0;
    for (int i = 0; i < n; i++) {
        list[i] = 0;
    }

    // get keywords containing polarization-related options
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        string = record.substr(next);
        iss.clear();
        iss.str(string);
        if (keyword == "POLARIZABLE") {
            int listInt;
            int j = nlist;
            while ((iss >> listInt) and (j < n)) {
                list[j] = listInt;
                j++;
            }
            while ((list[nlist] != 0) and (nlist <= n)) {
                nlist++;
            }
        }
        else if (keyword == "POLAR-PRINT") {
            polprt = true;
        }
        else if (keyword == "PCG-PRECOND") {
            pcgprec = true;
        }
        else if (keyword == "PCG-NOPRECOND") {
            pcgprec = false;
        }
        else if (keyword == "PCG-GUESS") {
            pcgguess = true;
        }
        else if (keyword == "PCG-NOGUESS") {
            pcgguess = false;
        }
        else if (keyword == "PCG-PEEK") {
            if (!(iss >> pcgpeek)) goto label_20;
        }
        else if (keyword == "TCG-GUESS") {
            tcgguess = true;
        }
        else if (keyword == "TCG-NOGUESS") {
            tcgguess = false;
        }
        else if (keyword == "TCG-PEEK") {
            if (!(iss >> tcgpeek)) goto label_20;
        }
        else if (keyword == "OPT-COEFF") {
            for (int j = 0; j <= maxopt; j++) {
                copt[j] = 0.;
            }
            real coptVal;
            int j = 0;
            while ((iss >> coptVal) and (j <= maxopt)) {
                copt[j] = coptVal;
                j++;
            }
        }
        label_20:;
    }

    // get maximum coefficient order for OPT induced dipoles
    if (poltyp.substr(0,3) == "OPT") {
        poltyp = "OPT";
        for (int i = 1; i <= maxopt; i++) {
            if (copt[i] != 0.) optorder = std::max(i,optorder);
        }
        for (int i = 0; i <= optorder; i++) {
            for (int j = optorder; j >= i; j--) {
               copm[i] += copt[j];
            }
        }
    }

    // perform dynamic allocation of some global arrays
    ipolar.allocate(n);
    polarity.allocate(n);
    thole.allocate(n);
    tholed.allocate(n);
    pdamp.allocate(n);
    udir.allocate(n);
    udirp.allocate(n);
    uind.allocate(n);
    uinp.allocate(n);
    douind.allocate(n);
    if (poltyp == "OPT") {
        uopt.allocate(n*(optorder+1));
        uoptp.allocate(n*(optorder+1));
        fopt.allocate(n*(optorder+1));
        foptp.allocate(n*(optorder+1));
    }

    // set the atoms allowed to have nonzero induced dipoles
    for (int i = 0; i < n; i++) {
        douind[i] = true;
    }
    i = 0;
    while ((list[i] != 0) and (i<n)) {
        if (i == 0) {
            for (int j = 0; j < n; j++) {
                douind[j] = false;
            }
        }
        if (list[i]>0 and list[i]<=n) {
            j = list[i] - 1;
            if (!douind[j]) {
                douind[j] = true;
            }
        }
        else if (list[i]<0 and list[i]>=-n) {
            for (int j = std::abs(list[i])-1; j < std::abs(list[i+1]); j++) {
                if (!douind[j]) {
                    douind[j] = true;
                }
            }
            i++;
        }
        i++;
    }

    // process keywords containing polarizability parameters
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "POLARIZE") {
            k = 0;
            pol = 0.;
            thl = -1.;
            thd = -1.;
            for (int j = 0; j < maxval; j++) {
                pg[j] = 0;
            }
            getnumb(record,k,next);
            gettext(record,text,next);
            iss.clear();
            iss.str(text);
            if (!(iss >> pol)) goto label_30;
            gettext(record,text,next);
            j = 0;
            getnumb(text,pg[0],j);
            if (pg[0] == 0) {
                iss.clear();
                iss.str(text);
                if (!(iss >> thl)) goto label_30;
                gettext(record,text,next);
                j = 0;
                getnumb(text,pg[0],j);
                string = record.substr(next);
                if (pg[0] == 0) {
                    iss.clear();
                    iss.str(text);
                    if (!(iss >> thd)) goto label_30;
                    iss.clear();
                    iss.str(string);
                    for (int j = 0; j < maxval; j++) {
                        if (!(iss >> pg[j])) goto label_30;
                    }
                }
                else {
                    iss.clear();
                    iss.str(string);
                    for (int j = 1; j < maxval; j++) {
                        if (!(iss >> pg[j])) goto label_30;
                    }
                }
            }
            else {
                string = record.substr(next);
                iss.clear();
                iss.str(string);
                for (int j = 1; j < maxval; j++) {
                    if (!(iss >> pg[j])) goto label_30;
                }
            }
            label_30:
            if (k > 0) {
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Atomic Dipole Polarizability Parameters :\n");
                    if (thd >= 0.) {
                        printf("\n     Atom Type           Alpha");
                        printf("       Thole      TholeD     Group Atom Types\n\n");
                    }
                    else if (thl >= 0.) {
                        printf("\n     Atom Type           Alpha");
                        printf("       Thole     Group Atom Types\n\n");
                    }
                    else {
                        printf("\n     Atom Type           Alpha");
                        printf("     Group Atom Types\n\n");
                    }
                }
                if (k <= maxtyp) {
                    int km = k-1;
                    polr[km] = pol;
                    athl[km] = std::max(0.,thl);
                    dthl[km] = std::max(0.,thd);
                    for (int j = 0; j < maxval; j++) {
                        pgrp[km][j] = pg[j]-1;
                        if (pg[j] == 0) {
                            npg = j;
                            break;
                        }
                    }
                    if (!silent) {
                        if (thd >= 0.) {
                            printf("    %8d        %10.3f  %10.3f  %10.3f       ", k,pol,thl,thd);
                            for (int j = 0; j < npg; j++) {
                                printf("%5d", pg[j]);
                            }
                            printf("\n");
                        }
                        else if (thl >= 0.) {
                            printf("    %8d        %10.3f  %10.3f       ", k,pol,thl);
                            for (int j = 0; j < npg; j++) {
                                printf("%5d", pg[j]);
                            }
                            printf("\n");
                        }
                        else {
                            printf("    %8d        %10.3f       ", k,pol);
                            for (int j = 0; j < npg; j++) {
                                printf("%5d", pg[j]);
                            }
                            printf("\n");
                        }
                    }
                }
                else {
                    printf("\n KPOLAR  --  Too many Dipole Polarizability Parameters\n");
                    informAbort = true;
                }
            }
        }
    }

    // process keywords with specific pair polarization values
    blank = "";
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "POLPAIR") {
            ia = 0;
            ib = 0;
            thl = -1.;
            thd = -1.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> thl >> thd;
            if (header and !silent) {
                header = false;
                printf("\n Additional Polarization Parameters for Specific Pairs :\n");
                if (thd >= 0.) {
                    printf("\n     Atom Types              Thole         TholeD\n\n");

                }
                else if (thl >= 0.) {
                    printf("\n     Atom Types              Thole\n\n");
                }
            }
            if (thd>=0. and !silent) {
                printf("      %4d%4d     %15.4f%15.4f\n", ia,ib,thl,thd);
            }
            else if (thl>=0. and !silent) {
                printf("      %4d%4d     %15.4f\n", ia,ib,thl);
            }
            size = 4;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
                pt = pa + pb;
            }
            else {
                pt = pb + pa;
            }
            for (int k = 0; k < maxnpp; k++) {
                if (kppr[k]==blank or kppr[k]==pt) {
                    kppr[k] = pt;
                    thlpr[k] = std::max(thl,0.);
                    thdpr[k] = std::max(thd,0.);
                    goto label_200;
                }
            }
            printf("\n KPOLAR  --  Too many Special Pair Thole Parameters\n");
            informAbort = true;
            label_200:;
        }
    }

    // find and store the atomic dipole polarizability parameters
    sixth = 1. / 6.;
    npolar = n;
    for (int i = 0; i < n; i++) {
        polarity[i] = 0.;
        thole[i] = 0.;
        tholed[i] = 0.;
        pdamp[i] = 0.;
        it = type[i];
        if (it != -1) {
            polarity[i] = polr[it];
            thole[i] = athl[it];
            tholed[i] = dthl[it];
            pdamp[i] = std::pow(polarity[i],sixth);
        }
    }

    // perform dynamic allocation of some global arrays
    jpolar.allocate(n);

    // perform dynamic allocation of some local arrays
    list.resize(n);
    rlist.resize(maxtyp);

    // set atom type index into condensed pair Thole matrices
    nlist = n;
    for (int i = 0; i < n; i++) {
        list[i] = type[i];
        jpolar[i] = list[i];
    }
    sortUnique(nlist, list, 0);
    for (int i = 0; i < maxtyp; i++) {
        rlist[i] = -1;
    }
    for (int i = 0; i < n; i++) {
        j = jpolar[i];
        if (rlist[j] == -1) {
            for (int k = 0; k < nlist; k++) {
                if (list[k] == j) rlist[j] = k;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        jpolar[i] = rlist[type[i]];
    }

    // perform dynamic allocation of some global arrays
    thlval.allocate(nlist*nlist);
    thdval.allocate(nlist*nlist);

    // use combination rules for pairwise Thole damping values
    for (int ii = 0; ii < nlist; ii++) {
        i = list[ii];
        for (int kk = ii; kk < nlist; kk++) {
            k = list[kk];
            thl = std::min(athl[i],athl[k]);
            if (thl == 0.) thl = std::max(athl[i],athl[k]);
            thd = std::min(dthl[i],dthl[k]);
            if (thd == 0.) thd = std::max(dthl[i],dthl[k]);
            thlval[kk*nlist+ii] = thl;
            thlval[ii*nlist+kk] = thl;
            thdval[kk*nlist+ii] = thd;
            thdval[ii*nlist+kk] = thd;
        }
    }

    // apply Thole damping values for special atom type pairs
    for (int i = 0; i < maxnpp; i++) {
        if (kppr[i] == blank) break;
        ia = rlist[std::stoi(kppr[i].substr(0,4))-1];
        ib = rlist[std::stoi(kppr[i].substr(4,4))-1];
        if (ia!=-1 and ib!=-1) {}
            thlval[ib*nlist+ia] = thlpr[i];
            thlval[ia*nlist+ib] = thlpr[i];
            thdval[ib*nlist+ia] = thdpr[i];
            thdval[ia*nlist+ib] = thdpr[i];
    }

    // setup exchange polarization via variable polarizability
    // kexpol(); // TODO

    // remove zero or undefined electrostatic sites from the list
    if ((use_polar or use_repel or use_solv) and !use_chgtrn) {
        npole = 0;
        ncp = 0;
        npolar = 0;
        nexpol = 0;
        for (int i = 0; i < n; i++) {
            if (polarity[i] == 0.) douind[i] = false;
            if (polsiz[i]!=0 or polarity[i]!=0.) {
                ipole[npole] = i;
                pollist[i] = npole;
                npole++;
                mono0[i] = pole[i][0];
                if (palpha[i] != 0.) ncp++;
                if (polarity[i] != 0.) {
                    ipolar[npolar] = i;
                    npolar++;
                    douind[i] = true;
                }
                if (tholed[i] != 0.) use_tholed = true;
                // if (kpep[i] != 0.) nexpol++;
            }
        }
    }

    // test multipoles at chiral sites and invert if necessary
    if (use_polar and !use_chgtrn) chkpole();

    // assign polarization group connectivity of each atom
    polargrp();

    // turn off polarizable multipole potentials if not used
    if (npole == 0) use_mpole = false;
    if (ncp != 0) use_chgpen = true;
    if (npolar == 0) use_polar = false;
    if (ncp != 0) use_thole = false;
    if (use_tholed) use_thole = true;
    if (nexpol != 0) use_expol = true;
}

/////////////////////////////////////////////////////
//                                                 //
//  polargrp  --  polarization group connectivity  //
//                                                 //
/////////////////////////////////////////////////////

// c     "polargrp" generates members of the polarization group of
// c     each atom and separate lists of the 1-2, 1-3 and 1-4 group
// c     connectivities

void polargrp()
{
    int i,j,k,m;
    int it,jt;
    int jj,kk;
    int start,stop;
    int nkeep,nlist;
    int maxkeep,maxlist;
    std::vector<int> keep;
    std::vector<int> list;
    std::vector<int> mask;
    bool done,qcmdAbort;

    // perform dynamic allocation of some global arrays
    np11.allocate(n);
    np12.allocate(n);
    np13.allocate(n);
    np14.allocate(n);
    ip11.allocate(n);
    ip12.allocate(n);
    ip13.allocate(n);
    ip14.allocate(n);

    // initialize size and connectivity of polarization groups
    for (int i = 0; i < n; i++) {
        np11[i] = 1;
        ip11[i][0] = i;
        np12[i] = 0;
        np13[i] = 0;
        np14[i] = 0;
    }

    // set termination flag and temporary group storage
    qcmdAbort = false;
    maxkeep = 100;
    maxlist = 10000;

    // find the directly connected group members for each atom
    for (int i = 0; i < n; i++) {
        it = type[i];
        if (it != -1) {
            for (int j = 0; j < n12[i]; j++) {
                jj = i12[i][j];
                jt = type[jj];
                for (int k = 0; k < maxval; k++) {
                    kk = pgrp[it][k];
                    if (kk == -1) goto label_20;
                    if (pgrp[it][k] == jt) {
                        if (np11[i] < maxp11) {
                            ip11[i][np11[i]] = jj;
                            np11[i]++;
                        }
                        else {
                            printf("\n POLARGRP  --  Too many Atoms in Polarization Group\n");
                            qcmdAbort = true;
                            goto label_30;
                        }
                    }
                }
                label_20:;
            }
        }
    }
    label_30:

    // make sure all connected group members are bidirectional
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < np11[i]; j++) {
            k = ip11[i][j];
            for (int m = 0; m < np11[k]; m++) {
                if (ip11[k][m] == i) goto label_50;
            }
            printf("\n POLARGRP  --  Check Polarization Groups for Atoms%9d and%9d\n", std::min(i,k)+1, std::max(i,k)+1);
            qcmdAbort = true;
            label_50:;
        }
    }

    // perform dynamic allocation of some local arrays
    keep.resize(maxkeep);
    list.resize(maxlist);
    mask.resize(n);

    // find any other group members for each atom in turn
    for (int i = 0; i < n; i++) {
        mask[i] = -1;
    }
    for (int i = 0; i < n; i++) {
        done = false;
        start = 0;
        stop = np11[i];
        for (int j = start; j < stop; j++) {
            jj = ip11[i][j];
            if (jj < i) {
                done = true;
                np11[i] = np11[jj];
                for (int k = 0; k < np11[i]; k++) {
                    ip11[i][k] = ip11[jj][k];
                }
            }
            else {
                mask[jj] = i;
            }
        }
        while (!done) {
            done = true;
            for (int j = start; j < stop; j++) {
                jj = ip11[i][j];
                for (int k = 0; k < np11[jj]; k++) {
                    kk = ip11[jj][k];
                    if (mask[kk] != i) {
                        if (np11[i] < maxp11) {
                            ip11[i][np11[i]] = kk;
                            np11[i]++;
                        }
                        else {
                            printf("\n POLARGRP  --  Too many Atoms in Polarization Group\n");
                            qcmdAbort = true;
                            goto label_70;
                        }
                        mask[kk] = i;
                    }
                }
            }
            if (np11[i] != stop) {
                done = false;
                start = stop;
                stop = np11[i];
            }
        }
        std::sort(ip11[i], ip11[i] + np11[i]);
    }
    label_70:
    if (qcmdAbort) fatal();

    // loop over atoms finding all the 1-2 group relationships
    for (int i = 0; i < n; i++) {
        mask[i] = -1;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < np11[i]; j++) {
            jj = ip11[i][j];
            mask[jj] = i;
        }
        nkeep = 0;
        for (int j = 0; j < np11[i]; j++) {
            jj = ip11[i][j];
            for (int k = 0; k < n12[jj]; k++) {
                kk = i12[jj][k];
                if (mask[kk] != i) {
                    keep[nkeep] = kk;
                    nkeep++;
                }
            }
        }
        nlist = 0;
        for (int j = 0; j < nkeep; j++) {
            jj = keep[j];
            for (int k = 0; k < np11[jj]; k++) {
                kk = ip11[jj][k];
                list[nlist] = kk;
                nlist++;
            }
        }
        sortUnique(nlist, list, 0);
        if (nlist <= maxp12) {
            np12[i] = nlist;
            for (int j = 0; j < nlist; j++) {
                ip12[i][j] = list[j];
            }
        }
        else {
            printf("\n POLARGRP  --  Too many Atoms in 1-2 Polarization Group\n");
            qcmdAbort = true;
            goto label_90;
        }
    }
    label_90:

    // loop over atoms finding all the 1-3 group relationships
    for (int i = 0; i < n; i++) {
        mask[i] = -1;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < np11[i]; j++) {
            jj = ip11[i][j];
            mask[jj] = i;
        }
        for (int j = 0; j < np12[i]; j++) {
            jj = ip12[i][j];
            mask[jj] = i;
        }
        nlist = 0;
        for (int j = 0; j < np12[i]; j++) {
            jj = ip12[i][j];
            for (int k = 0; k < np12[jj]; k++) {
                kk = ip12[jj][k];
                if (mask[kk] != i) {
                    list[nlist] = kk;
                    nlist++;
                }
            }
        }
        sortUnique(nlist, list, 0);
        if (nlist <= maxp13) {
            np13[i] = nlist;
            for (int j = 0; j < nlist; j++) {
                ip13[i][j] = list[j];
            }
        }
        else {
            printf("\n POLARGRP  --  Too many Atoms in 1-3 Polarization Group\n");
            qcmdAbort = true;
            goto label_110;
        }
    }
    label_110:

    // loop over atoms finding all the 1-4 group relationships
    for (int i = 0; i < n; i++) {
        mask[i] = -1;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < np11[i]; j++) {
            jj = ip11[i][j];
            mask[jj] = i;
        }
        for (int j = 0; j < np12[i]; j++) {
            jj = ip12[i][j];
            mask[jj] = i;
        }
        for (int j = 0; j < np13[i]; j++) {
            jj = ip13[i][j];
            mask[jj] = i;
        }
        nlist = 0;
        for (int j = 0; j < np13[i]; j++) {
            jj = ip13[i][j];
            for (int k = 0; k < np12[jj]; k++) {
                kk = ip12[jj][k];
                if (mask[kk] != i) {
                    list[nlist] = kk;
                    nlist++;
                }
            }
        }
        sortUnique(nlist, list, 0);
        if (nlist <= maxp14) {
            np14[i] = nlist;
            for (int j = 0; j < nlist; j++) {
                ip14[i][j] = list[j];
            }
        }
        else {
            printf("\n POLARGRP  --  Too many Atoms in 1-4 Polarization Group\n");
            qcmdAbort = true;
            goto label_130;
        }
    }
    label_130:
    if (qcmdAbort) fatal();
}
}
