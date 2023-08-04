//////////////////////////////////////////////////////
//                                                  //
//  kmpole.cpp  --  multipole parameter assignment  //
//                                                  //
//////////////////////////////////////////////////////

// "kmpole" assigns atomic multipole moments to the atoms of
// the structure and processes any new or changed values


#include "atomid.h"
#include "atoms.h"
#include "chgpen.h"
#include "chkpole.h"
#include "couple.h"
#include "gettext.h"
#include "inform.h"
#include "kcpen.h"
#include "keys.h"
#include "kmpole.h"
#include "kmulti.h"
#include "mathConst.h"
#include "mplpot.h"
#include "mpole.h"
#include "numeral.h"
#include "polar.h"
#include "polgrp.h"
#include "potent.h"
#include "units.h"
#include "upcase.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

void kmpole()
{
    int i,j,k,l,m;
    int ji,ki,li;
    int it,jt,kt,lt;
    int ic,imp,nmp;
    int size,next;
    int number;
    int kz,kx,ky;
    int ztyp,xtyp,ytyp;
    int polmax;
    std::vector<int> mpt;
    std::vector<int> mpz;
    std::vector<int> mpx;
    std::vector<int> mpy;
    double pel,pal;
    double mpl[13];
    bool header,path;
    std::string pa,pb,pc,pd;
    std::string axt;
    std::string blank,pt;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // count the number of existing multipole parameters
    blank = "";
    nmp = maxnmp;
    for (int i = maxnmp-1; i >=0; i--) {
        if (kmp[i] == blank)  nmp = i;
    }

    // find and count new multipole parameters in the keyfile
    imp = 0;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "MULTIPOLE") {
            k = 0;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            int numObjects = 0;
            int value;
            while (iss >> value) {
                numObjects++;
            }
            iss.clear();
            iss.seekg(0);
            if (numObjects == 5) {
                if (!(iss >> k >> kz >> kx >> ky >> mpl[0])) goto label_50;
            }
            else if (numObjects == 4) {
                if (!(iss >> k >> kz >> kx >> mpl[0])) goto label_50;
            }
            else if (numObjects == 3) {
                if (!(iss >> k >> kz >> mpl[0])) goto label_50;
            }
            else if (numObjects == 2) {
                if (!(iss >> k >> mpl[0])) goto label_50;
            }
            else {
                goto label_50;
            }
            if (k > 0) {
                record = keyline[i+1];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[1] >> mpl[2] >> mpl[3])) goto label_50;
                record = keyline[i+2];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[4])) goto label_50;
                record = keyline[i+3];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[7] >> mpl[8])) goto label_50;
                record = keyline[i+4];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[10] >> mpl[11] >> mpl[12])) goto label_50;
                imp++;
            }
        label_50:
        continue;
        }
    }

    // check for too many combined parameter values
    nmp += imp;
    if (nmp > maxnmp) {
        printf("\n KMPOLE  --  Too many Atomic Multipole Parameters\n");
        informAbort = true;
    }

    // move existing parameters to make room for new values
    if (imp != 0) {
        for (int j = nmp-1; j >= imp; j--) {
            k = j - imp;
            kmp[j] = kmp[k];
            mpaxis[j] = mpaxis[k];
            for (int m = 0; m < 13; m++) {
                multip[j][m] = multip[k][m];
            }
        }
    }

    // process keywords containing atomic multipole parameters
    imp = 0;
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "MULTIPOLE") {
            k = 0;
            kz = 0;
            kx = 0;
            ky = 0;
            axt = "Z-then-X";
            for (int j = 0; j < 13; j++) {
                mpl[j] = 0.;
            }
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            int numObjects = 0;
            int value;
            while (iss >> value) {
                numObjects++;
            }
            iss.clear();
            iss.seekg(0);
            if (numObjects == 5) {
                if (!(iss >> k >> kz >> kx >> ky >> mpl[0])) goto label_130;
            }
            else if (numObjects == 4) {
                if (!(iss >> k >> kz >> kx >> mpl[0])) goto label_130;
            }
            else if (numObjects == 3) {
                if (!(iss >> k >> kz >> mpl[0])) goto label_130;
            }
            else if (numObjects == 2) {
                if (!(iss >> k >> mpl[0])) goto label_130;
            }
            else {
                goto label_130;
            }
            if (k > 0) {
                if (kz == 0) axt = "None";
                if (kz!=0 and kx==0) axt = "Z-Only";
                if (kz<0 or kx<0) axt = "Bisector";
                if (kx<0 and ky<0) axt = "Z-Bisect";
                if (std::max({kz,kx,ky}) < 0) axt = "3-Fold";
                kz = std::abs(kz);
                kx = std::abs(kx);
                ky = std::abs(ky);
                record = keyline[i+1];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[1] >> mpl[2] >> mpl[3])) goto label_130;
                record = keyline[i+2];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[4])) goto label_130;
                record = keyline[i+3];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[7] >> mpl[8])) goto label_130;
                record = keyline[i+4];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[10] >> mpl[11] >> mpl[12])) goto label_130;
                mpl[5] = mpl[7];
                mpl[6] = mpl[10];
                mpl[9] = mpl[11];
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Atomic Multipole Parameters :\n\n");
                    printf("     Atom Type     Coordinate Frame");
                    printf(" Definition         Multipole Moments\n");
                }
                if (!silent) {
                    printf("\n      %6d   %6d %6d %6d   %-8s   %9.5f\n", k, kz, kx, ky, axt.c_str(), mpl[0]);
                    printf("                                                 %9.5f%9.5f%9.5f\n", mpl[1], mpl[2], mpl[3]);
                    printf("                                                 %9.5f\n", mpl[4]);
                    printf("                                                 %9.5f%9.5f\n", mpl[7], mpl[8]);
                    printf("                                                 %9.5f%9.5f%9.5f\n", mpl[10], mpl[11], mpl[12]);
                }
                size = 4;
                pa = numeral (k,size);
                pb = numeral (kz,size);
                pc = numeral (kx,size);
                pd = numeral (ky,size);
                pt = pa + pb + pc + pd;
                kmp[imp] = pt;
                mpaxis[imp] = axt;
                for (int j = 0; j < 13; j++) {
                    multip[imp][j] = mpl[j];
                }
                imp++;
            }
            label_130:
            continue;
        }
    }

    // perform dynamic allocation of some global arrays
    if (ipole.size() != 0) ipole.resize(0);
    if (polsiz.size() != 0) polsiz.resize(0);
    if (pollist.size() != 0) pollist.resize(0);
    if (zaxis.size() != 0) zaxis.resize(0);
    if (xaxis.size() != 0) xaxis.resize(0);
    if (yaxis.size() != 0) yaxis.resize(0);
    if (pole.size() != 0) pole.resize(0);
    if (rpole.size() != 0) rpole.resize(0);
    if (mono0.size() != 0) mono0.resize(0);
    if (polaxe.size() != 0) polaxe.resize(0);
    if (np11.size() != 0) np11.resize(0);
    if (np12.size() != 0) np12.resize(0);
    if (np13.size() != 0) np13.resize(0);
    if (np14.size() != 0) np14.resize(0);
    ipole.resize(n, -1);
    polsiz.resize(n, 0);
    pollist.resize(n, -1);
    zaxis.resize(n, 0);
    xaxis.resize(n, 0);
    yaxis.resize(n, 0);
    pole.resize(n, std::vector<double>(maxpole, 0.));
    rpole.resize(n, std::vector<double>(maxpole));
    mono0.resize(n, 0.);
    polaxe.resize(n, "None");
    np11.resize(n, 0);
    np12.resize(n, 0);
    np13.resize(n, 0);
    np14.resize(n, 0);

    // perform dynamic allocation of some local arrays
    mpt.resize(maxnmp);
    mpz.resize(maxnmp);
    mpx.resize(maxnmp);
    mpy.resize(maxnmp);

    // store the atom types associated with each parameter
    for (int i = 0; i < nmp; i++) {
        mpt[i] = std::stoi(kmp[i].substr(0,4)) - 1;
        mpz[i] = std::stoi(kmp[i].substr(4,4)) - 1;
        mpx[i] = std::stoi(kmp[i].substr(8,4)) - 1;
        mpy[i] = std::stoi(kmp[i].substr(12,4)) - 1;
    }

    // assign multipole parameters via only 1-2 connected atoms
    for (int i = 0; i < n; i++) {
        it = type[i];
        for (int imp = 0; imp < nmp; imp++) {
            if (it == mpt[imp]) {
                ztyp = mpz[imp];
                xtyp = mpx[imp];
                ytyp = mpy[imp];
                for (int j = 0; j < n12[i]; j++) {
                    ji = i12[i][j];
                    jt = type[ji];
                    if (jt == ztyp) {
                        for (int k = 0; k < n12[i]; k++) {
                            ki = i12[i][k];
                            kt = type[ki];
                            if (kt==xtyp and ki!=ji) {
                                if (ytyp == -1) {
                                    pollist[i] = i;
                                    zaxis[i] = ji+1;
                                    xaxis[i] = ki+1;
                                    polaxe[i] = mpaxis[imp];
                                    for (int m = 0; m < 13; m++) {
                                        pole[i][m] = multip[imp][m];
                                    }
                                    goto label_140;
                                }
                                for (int l = 0; l < n12[i]; l++) {
                                    li = i12[i][l];
                                    lt = type[li];
                                    if (lt==ytyp and li!=ji and li!=ki) {
                                        pollist[i] = i;
                                        zaxis[i] = ji+1;
                                        xaxis[i] = ki+1;
                                        yaxis[i] = li+1;
                                        polaxe[i] = mpaxis[imp];
                                        for (int m = 0; m < 13; m++) {
                                            pole[i][m] = multip[imp][m];
                                        }
                                        goto label_140;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // assign multipole parameters via 1-2 and 1-3 connected atoms
        for (int imp = 0; imp < nmp; imp++) {
            if (it == mpt[imp]) {
                ztyp = mpz[imp];
                xtyp = mpx[imp];
                ytyp = mpy[imp];
                for (int j = 0; j < n12[i]; j++) {
                    ji = i12[i][j];
                    jt = type[ji];
                    if (jt == ztyp) {
                        for (int k = 0; k < n13[i]; k++) {
                            ki = i13[i][k];
                            kt = type[ki];
                            path = false;
                            for (int m = 0; m < n12[ki]; m++) {
                                if (i12[ki][m] == ji) path = true;
                            }
                            if (kt==xtyp and path) {
                                if (ytyp == -1) {
                                    pollist[i] = i;
                                    zaxis[i] = ji+1;
                                    xaxis[i] = ki+1;
                                    polaxe[i] = mpaxis[imp];
                                    for (int m = 0; m < 13; m++) {
                                        pole[i][m] = multip[imp][m];
                                    }
                                    goto label_140;
                                }
                                for (int l = 0; l < n13[i]; l++) {
                                    li = i13[i][l];
                                    lt = type[li];
                                    path = false;
                                    for (int m = 0; m < n12[li]; m++) {
                                        if (i12[li][m] == ji) path = true;
                                    }
                                    if (lt==ytyp and li!=ki and path) {
                                        pollist[i] = i;
                                        zaxis[i] = ji+1;
                                        xaxis[i] = ki+1;
                                        yaxis[i] = li+1;
                                        polaxe[i] = mpaxis[imp];
                                        for (int m = 0; m < 13; m++) {
                                            pole[i][m] = multip[imp][m];
                                        }
                                        goto label_140;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // assign multipole parameters via only a z-defining atom
        for (int imp = 0; imp < nmp; imp++) {
            if (it == mpt[imp]) {
                ztyp = mpz[imp];
                xtyp = mpx[imp];
                ytyp = mpy[imp];
                for (int j = 0; j < n12[i]; j++) {
                    ji = i12[i][j];
                    jt = type[ji];
                    if (jt == ztyp) {
                        if (xtyp == -1) {
                            pollist[i] = i;
                            zaxis[i] = ji+1;
                            polaxe[i] = mpaxis[imp];
                            for (int m = 0; m < 13; m++) {
                                pole[i][m] = multip[imp][m];
                            }
                            goto label_140;
                        }
                    }
                }
            }
        }

        // assign multipole parameters via no connected atoms
        for (int imp = 0; imp < nmp; imp++) {
            if (it == mpt[imp]) {
                ztyp = mpz[imp];
                xtyp = mpx[imp];
                ytyp = mpy[imp];
                if (ztyp == -1) {
                    pollist[i] = i;
                    polaxe[i] = mpaxis[imp];
                    for (int m = 0; m < 13; m++) {
                        pole[i][m] = multip[imp][m];
                    }
                    goto label_140;
                }
            }
        }
        label_140:
        continue;
    }

    // process keywords with multipole parameters for specific atoms
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "MULTIPOLE") {
            k = 0;
            kz = 0;
            kx = 0;
            ky = 0;
            axt = "Z-then-X";
            for (int j = 0; j < 13; j++) {
                mpl[j] = 0.;
            }
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            int numObjects = 0;
            int value;
            while (iss >> value) {
                numObjects++;
            }
            iss.clear();
            iss.seekg(0);
            if (numObjects == 5) {
                if (!(iss >> k >> kz >> kx >> ky >> mpl[0])) goto label_210;
            }
            else if (numObjects == 4) {
                if (!(iss >> k >> kz >> kx >> mpl[0])) goto label_210;
            }
            else if (numObjects == 3) {
                if (!(iss >> k >> kz >> mpl[0])) goto label_210;
            }
            else if (numObjects == 2) {
                if (!(iss >> k >> mpl[0])) goto label_210;
            }
            else {
                goto label_210;
            }
            if (k<0 and k>=-n) {
                k = -k;
                if (kz == 0)  axt = "None";
                if (kz!=0 and kx==0) axt = "Z-Only";
                if (kz<0 or kx<0) axt = "Bisector";
                if (kx<0 and ky<0) axt = "Z-Bisect";
                if (std::max({kz,kx,ky}) < 0) axt = "3-Fold";
                kz = std::abs(kz);
                kx = std::abs(kx);
                ky = std::abs(ky);
                record = keyline[i+1];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[1] >> mpl[2] >> mpl[3])) goto label_210;
                record = keyline[i+2];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[4])) goto label_210;
                record = keyline[i+3];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[7] >> mpl[8])) goto label_210;
                record = keyline[i+4];
                iss.clear();
                iss.str(record);
                if (!(iss >> mpl[10] >> mpl[11] >> mpl[12])) goto label_210;
                mpl[5] = mpl[7];
                mpl[6] = mpl[10];
                mpl[9] = mpl[11];
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Atomic Multipoles for Specific Atoms :\n\n");
                    printf("     Atom          Coordinate Frame");
                    printf(" Definition         Multipole Moments\n");
                }
                if (!silent) {
                    printf("\n   %6d      %6d %6d %6d   %-8s   %9.5f\n", k, kz, kx, ky, axt.c_str(), mpl[0]);
                    printf("                                                 %9.5f%9.5f%9.5f\n", mpl[1], mpl[2], mpl[3]);
                    printf("                                                 %9.5f\n", mpl[4]);
                    printf("                                                 %9.5f%9.5f\n", mpl[7], mpl[8]);
                    printf("                                                 %9.5f%9.5f%9.5f\n", mpl[10], mpl[11], mpl[12]);
                }
                k--;
                kz--;
                kx--;
                ky--;
                pollist[k] = k;
                zaxis[k] = kz+1;
                xaxis[k] = kx+1;
                yaxis[k] = ky+1;
                polaxe[k] = axt;
                for (int j = 0; j < 13; j++) {
                    pole[k][j] = mpl[j];
                }
            }
            label_210:
            continue;
        }
    }

    // convert the dipole and quadrupole moments to Angstroms,
    // quadrupole divided by 3 for use as traceless values
    for (int i = 0; i < n; i++) {
        for (int k = 1; k < 4; k++) {
            pole[i][k] = pole[i][k] * bohr;
        }
        for (int k = 4; k < 13; k++) {
            pole[i][k] = pole[i][k] * bohr*bohr / 3.;
        }
    }

    // get the order of the multipole expansion at each site
    npole = n;
    polmax = 0;
    for (int i = 0; i < n; i++) {
        size = 0;
        for (int k = 0; k < maxpole; k++) {
            if (pole[i][k] != 0.)  size = std::max(k+1,size);
        }
        if (size > 4) size = 13;
        else if (size > 1) size = 4;
        polsiz[i] = size;
        polmax = std::max(polmax,size);
    }

    // warn if there are sites with no atomic multipole values
    if (polmax != 0) {
        header = true;
        for (int i = 0; i < n; i++) {
            if (pollist[i] == -1) {
                if (header) {
                    header = false;
                    printf("\n Undefined Atomic Multipole Parameters :\n\n");
                }
                printf(" Warning, No Multipole Parameters for Atom%7d\n", i+1);
            }
            pollist[i] = -1;
        }
    }

    // perform dynamic allocation of some global arrays
    // if polarization not used, zero out induced dipoles
    if (!use_polar) {
        if (uind.size() != 0) uind.resize(0);
        if (uinp.size() != 0) uinp.resize(0);
        if (uinds.size() != 0) uinds.resize(0);
        if (uinps.size() != 0) uinps.resize(0);
        uind.resize(n, std::vector<double>(3, 0.));
        uinp.resize(n, std::vector<double>(3, 0.));
        uinds.resize(n, std::vector<double>(3, 0.));
        uinps.resize(n, std::vector<double>(3, 0.));
    }

    // perform dynamic allocation of some global arrays
    if (pcore.size() != 0) pcore.resize(0);
    if (pval.size() != 0) pval.resize(0);
    if (pval0.size() != 0) pval0.resize(0);
    if (palpha.size() != 0) palpha.resize(0);
    pcore.resize(n);
    pval.resize(n);
    pval0.resize(n);
    palpha.resize(n);

    // test
    // find new charge penetration parameters in the keyfile
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "CHGPEN") {
            k = 0;
            pel = 0.;
            pal = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            if (!(iss >> k >> pel >> pal)) goto label_260;
            if (k > 0) {
                k--;
                cpele[k] = std::abs(pel);
                cpalp[k] = pal;
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Charge Penetration Parameters :\n\n");
                    printf("     Atom Class           Core Chg           Damp\n\n");
                }
                if (!silent) {
                    printf("      %6d       %15.3f%15.4f\n", k+1, pel, pal);
                }
            }
            label_260:
            continue;
        }
    }

    // assign the charge penetration charge and alpha parameters 
    ncp = 0;
    for (int i = 0; i < n; i++) {
        pcore[i] = 0.;
        pval[i] = pole[i][0];
        pval0[i] = pval[i];
        palpha[i] = 0.;
        ic = atomClass[i];
        if (ic != 0) {
            pcore[i] = cpele[ic];
            pval[i] = pole[i][0] - cpele[ic];
            pval0[i] = pval[i];
            palpha[i] = cpalp[ic];
        }
    }

    // test
    // process keywords with charge penetration for specific atoms
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "CHGPEN") {
            k = 0;
            pel = 0.;
            pal = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            if (!(iss >> k >> pel >> pal)) goto label_290;
            if (k<0 and k>=-n) {
                k = -k;
                k--;
                pcore[k] = std::abs(pel);
                pval[k] = pole[k][0] - std::abs(pel);
                palpha[k] = pal;
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Charge Penetration for Specific Atoms :\n\n");
                    printf("     Atom                 Core Chg           Damp\n\n");
                }
                if (!silent) {
                    printf("      %6d       %15.3f%15.4f\n", k+1, pel, pal);
                }
            }
            label_290:
            continue;
        }
    }

    // remove zero or undefined electrostatic sites from the list
    if ((use_mpole or use_repel or use_solv) and !use_polar and !use_chgtrn) {
        npole = 0;
        ncp = 0;
        for (int i = 0; i < n; i++) {
            if (polsiz[i] != 0) {
                ipole[npole] = i;
                pollist[i] = npole;
                npole++;
                mono0[i] = pole[i][0];
                if (palpha[i] != 0.)  ncp++;
            }
        }
    }

    // test multipoles at chiral sites and invert if necessary
    if (use_mpole and !use_polar and !use_chgtrn) chkpole();

    // turn off atomic multipole potentials if not used
    if (npole == 0)  use_mpole = false;
    if (ncp != 0)  use_chgpen = true;
}
