// Author: Moses KJ Chung
// Year:   2023

#include "atomid.h"
#include "atoms.h"
#include "couple.h"
#include "fields.h"
#include "getnumb.h"
#include "gettext.h"
#include "inform.h"
#include "keys.h"
#include "khbond.h"
#include "kvdw.h"
#include "kvdws.h"
#include "kvdwpr.h"
#include "mathConst.h"
#include "merck.h"
#include "numeral.h"
#include "potent.h"
#include "sort.h"
#include "upcase.h"
#include "vdw.h"
#include "vdwpot.h"
#include <cmath>
#include <sstream>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  kvdw  --  van der Waals parameter assignment  //
//                                                //
////////////////////////////////////////////////////

// "kvdw" assigns the parameters to be used in computing the
// van der Waals interactions and processes any new or changed
// values for these parameters

void kvdw()
{
    int i,j,k;
    int ii,kk;
    int ia,ib;
    int next,size;
    int maxdim;
    int nlist,number;
    std::vector<int> list;
    real rd,ep,rdn,gik;
    std::vector<real> srad;
    std::vector<real> srad4;
    std::vector<real> seps;
    std::vector<real> seps4;
    bool header;
    std::string pa,pb;
    std::string blank,pt;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // process keywords containing van der Waals parameters
    maxdim = maxclass;
    if (vdwindex == "TYPE") maxdim = maxtyp;
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "VDW") {
            getnumb(record,k,next);
            if (k>0 and k<=maxdim) {
                int km = k-1;
                rd = 0.;
                ep = 0.;
                rdn = 0.;
                string = record.substr(next);
                iss.clear();
                iss.str(string);
                iss >> rd >> ep >> rdn;
                if (header and !silent) {
                    header = false;
                    if (vdwindex == "CLASS") {
                        printf("\n Additional van der Waals Parameters :\n\n");
                        printf("     Atom Class               Size        Epsilon        Reduction\n\n");
                    }
                    else {
                        printf("\n Additional van der Waals Parameters :\n\n");
                        printf("     Atom Type                Size        Epsilon        Reduction\n\n");
                    }
                }
                rad[km] = rd;
                eps[km] = ep;
                reduct[km] = rdn;
                if (!silent) {
                    printf("      %6d       %15.4f%15.4f%15.3f\n", k,rd,ep,rdn);
                }
            }
            else if (k > maxclass) {
                printf("\n KVDW  --  Only Atom Classes through%4d are Allowed\n", maxclass);
                informAbort = true;
            }
        }
    }

    // process keywords containing 1-4 van der Waals parameters
    maxdim = maxclass;
    if (vdwindex == "TYPE") maxdim = maxtyp;
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "VDW14") {
            getnumb(record,k,next);
            if (k>0 and k<=maxdim) {
                int km = k - 1;
                rd = 0.;
                ep = 0.;
                string = record.substr(next);
                iss.clear();
                iss.str(string);
                iss >> rd >> ep;
                if (header and !silent) {
                    header = false;
                    if (vdwindex == "CLASS") {
                        printf("\n Additional 1-4 van der Waals Parameters :\n\n");
                        printf("     Atom Class               Size        Epsilon\n\n");
                    }
                    else {
                        printf("\n Additional 1-4 van der Waals Parameters :\n\n");
                        printf("     Atom Type                Size        Epsilon\n\n");
                    }
                }
                rad4[km] = rd;
                eps4[km] = ep;
                if (!silent) {
                    printf("      %6d       %15.4f%15.4f\n", k, rd, ep);
                }
            }
            else if (k > maxclass) {
                printf("\n KVDW  --  Only Atom Classes through%4d are Allowed\n", maxclass);
               informAbort = true;
            }
        }
    }

    // process keywords containing specific pair vdw parameters
    blank = "";
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "VDWPAIR" or keyword == "VDWPR") {
            ia = 0;
            ib = 0;
            rd = 0.;
            ep = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            if (!(iss >> ia >> ib >> rd >> ep)) goto label_150;
            if (header and !silent) {
                header = false;
                if (vdwindex == "CLASS") {
                    printf("\n Additional van der Waals Parameters for Specific Pairs :\n\n");
                    printf("     Atom Classes         Size Sum        Epsilon\n\n");
                }
                else {
                    printf("\n Additional van der Waals Parameters for Specific Pairs :\n\n");
                    printf("     Atom Types           Size Sum        Epsilon\n\n");
                }
            }
            if (!silent) {
                printf("      %4d%4d     %15.4f%15.4f\n", ia, ib, rd, ep);
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
            for (int k = 0; k < maxnvp; k++) {
                if (kvpr[k]==blank or kvpr[k]==pt) {
                    kvpr[k] = pt;
                    radpr[k] = rd;
                    epspr[k] = ep;
                    goto label_150;
                }
            }
            printf("\n KVDW  --  Too many Special Pair VDW Parameters\n");
            informAbort = true;
            label_150:
            continue;
        }
    }

    // process keywords containing hydrogen bonding vdw parameters
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "HBOND") {
            ia = 0;
            ib = 0;
            rd = 0.;
            ep = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            if (!(iss >> ia >> ib >> rd >> ep)) goto label_200;
            if (header and !silent) {
                header = false;
                if (vdwindex == "CLASS") {
                    printf("\n Additional van der Waals Hydrogen Bonding Parameters :\n\n");
                    printf("     Atom Classes         Size Sum        Epsilon\n\n");
                    
                }
                else {
                    printf("\n Additional van der Waals Hydrogen Bonding Parameters :\n\n");
                    printf("     Atom Types           Size Sum        Epsilon\n\n");
                }
            }
            if (!silent) {
                printf("      %4d%4d     %15.4f%15.4f\n", ia, ib, rd, ep);
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
            for (int k = 0; k < maxnvp; k++) {
                if (khb[k]==blank or khb[k]==pt) {
                    khb[k] = pt;
                    radhb[k] = rd;
                    epshb[k] = ep;
                    goto label_200;
                }
            }
            printf("\n KVDW  --  Too many Hydrogen Bonding Pair VDW Parameters\n");
            informAbort = true;
            label_200:
            continue;
        }
    }

    // perform dynamic allocation of some global arrays
    if (ivdw.size() != 0) ivdw.resize(0);
    if (jvdw.size() != 0) jvdw.resize(0);
    if (mvdw.size() != 0) mvdw.resize(0);
    if (ired.size() != 0) ired.resize(0);
    if (kred.size() != 0) kred.resize(0);
    if (xred.size() != 0) xred.resize(0);
    if (yred.size() != 0) yred.resize(0);
    if (zred.size() != 0) zred.resize(0);
    ivdw.resize(n);
    jvdw.resize(n);
    mvdw.resize(maxtyp);
    ired.resize(n);
    kred.resize(n);
    xred.resize(n);
    yred.resize(n);
    zred.resize(n);

    // perform dynamic allocation of some local arrays
    list.resize(n);
    srad.resize(maxtyp);
    srad4.resize(maxtyp);
    seps.resize(maxtyp);
    seps4.resize(maxtyp);

    // set type or class index into condensed pair matrices
    nlist = n;
    for (int i = 0; i < n; i++) {
        list[i] = -1;
        if (vdwindex == "TYPE") {
            list[i] = type[i];
        }
        else {
            list[i] = atomClass[i];
        }
        jvdw[i] = list[i];
    }
    sortUnique(nlist,list,0);
    for (int i = 0; i < maxtyp; i++) {
        mvdw[i] = -1;
    }
    for (int i = 0; i < n; i++) {
        j = jvdw[i];
        if (mvdw[j] == -1) {
            for (int k = 0; k < nlist; k++) {
                if (list[k] == j) mvdw[j] = k;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        if (vdwindex == "TYPE") {
            k = type[i];
            jvdw[i] = mvdw[k];
        }
        else {
            k = atomClass[i];
            jvdw[i] = mvdw[k];
        }
    }

    // get the vdw radii and well depths for each atom type
    maxdim = maxclass;
    if (vdwindex == "TYPE") maxdim = maxtyp;
    for (int i = 0; i < maxdim; i++) {
        if (rad4[i] == 0.) rad4[i] = rad[i];
        if (eps4[i] == 0.) eps4[i] = eps[i];
        if (radtyp == "SIGMA") {
            rad[i] = twosix * rad[i];
            rad4[i] = twosix * rad4[i];
        }
        if (radsiz == "DIAMETER") {
            rad[i] = 0.5 * rad[i];
            rad4[i] = 0.5 * rad4[i];
        }
        srad[i] = std::sqrt(rad[i]);
        eps[i] = std::abs(eps[i]);
        seps[i] = std::sqrt(eps[i]);
        srad4[i] = std::sqrt(rad4[i]);
        eps4[i] = std::abs(eps4[i]);
        seps4[i] = std::sqrt(eps4[i]);
    }

    // perform dynamic allocation of some global arrays
    if (radmin.size() != 0) radmin.resize(0);
    if (epsilon.size() != 0) epsilon.resize(0);
    if (radmin4.size() != 0) radmin4.resize(0);
    if (epsilon4.size() != 0) epsilon4.resize(0);
    if (radhbnd.size() != 0) radhbnd.resize(0);
    if (epshbnd.size() != 0) epshbnd.resize(0);
    radmin.resize(nlist, std::vector<real>(nlist));
    epsilon.resize(nlist, std::vector<real>(nlist));
    radmin4.resize(nlist, std::vector<real>(nlist));
    epsilon4.resize(nlist, std::vector<real>(nlist));
    radhbnd.resize(nlist, std::vector<real>(nlist));
    epshbnd.resize(nlist, std::vector<real>(nlist));

    // use combination rules to set pairwise vdw radii sums
    for (int ii = 0; ii < nlist; ii++) {
        i = list[ii];
        for (int kk = ii; kk < nlist; kk++) {
            k = list[kk];
            if (radrule == "MMFF94") {
                if (i != k) {
                    rd = 0.5 * (rad[i]+rad[k]);
                    if (da[i]!='D' and da[k]!='D') {
                        if (rd != 0.) {
                            gik = (rad[i]-rad[k])/(rad[i]+rad[k]);
                            rd = (1.+0.2*(1.-exp(-12.*gik*gik))) * rd;
                        }
                    }
                }
                else rd = rad[i];
            }
            else if (rad[i]==0. and rad[k]==0.) rd = 0.;
            else if (radrule == "ARITHMETIC") rd = rad[i] + rad[k];
            else if (radrule == "GEOMETRIC") rd = 2. * (srad[i] * srad[k]);
            else if (radrule == "CUBIC-MEAN") {
                real radi = rad[i];
                real radk = rad[k];
                real radi2 = radi * radi;
                real radk2 = radk * radk;
                real radi3 = radi2 * radi;
                real radk3 = radk2 * radk;
                rd = 2. * (radi3+radk3)/(radi2+radk2);
            }
            else rd = rad[i] + rad[k];
            
            radmin[kk][ii] = rd;
            radmin[ii][kk] = rd;
        }
    }

    // use combination rules to set pairwise well depths
    for (int ii = 0; ii < nlist; ii++) {
        i = list[ii];
        for (int kk = ii; kk < nlist; kk++) {
            k = list[kk];
            if (epsrule == "MMFF94") {
                ep = 0.;
                if (nn[i]!=0. and nn[k]!=0. and radmin[kk][ii]!=0.) {
                    ep = 181.16*g[i]*g[k]*alph[i]*alph[k] / ((std::sqrt(alph[i]/nn[i])+std::sqrt(alph[k]/nn[k]))*std::pow(radmin[kk][ii],6));
                }
                if (i == k) eps[i] = ep;
            }
            else if (eps[i]==0. and eps[k]==0.) ep = 0.;
            else if (epsrule == "ARITHMETIC") ep = 0.5 * (eps[i] + eps[k]);
            else if (epsrule == "GEOMETRIC") ep = seps[i] * seps[k];
            else if (epsrule == "HARMONIC") ep = 2. * (eps[i]*eps[k]) / (eps[i]+eps[k]);
            else if (epsrule == "HHG") {
                real sepsik = seps[i]+seps[k];
                real sepsik2 = sepsik * sepsik;
                ep = 4. * (eps[i]*eps[k]) / sepsik2;
            }
            else if (epsrule == "W-H") {
                real radik = rad[i]*rad[k];
                real radik3 = std::pow(radik, 3);
                real radik6 = std::pow(rad[i],6) + std::pow(rad[k],6);
                ep = 2. * (seps[i]*seps[k]) * radik3 / radik6;
            }
            else ep = seps[i] * seps[k];
            epsilon[kk][ii] = ep;
            epsilon[ii][kk] = ep;
        }
    }

    // use combination rules to set pairwise 1-4 vdw radii sums
    for (int ii = 0; ii < nlist; ii++) {
        i = list[ii];
        for (int kk = ii; kk < nlist; kk++) {
            k = list[kk];
            if (radrule == "MMFF94") {
                if (i != k) {
                    rd = 0.5 * (rad[i]+rad[k]);
                    if (da[i]!='D' and da[k]!='D') {
                        if (rd != 0.) {
                            gik = (rad[i]-rad[k])/(rad[i]+rad[k]);
                            rd = (1.+0.2*(1.-std::exp(-12.*gik*gik))) * rd;
                        }
                    }
                }
                else rd = rad[i];
            }
            else if (rad4[i]==0. and rad4[k]==0.) rd = 0.;
            else if (radrule == "ARITHMETIC") rd = rad4[i] + rad4[k];
            else if (radrule == "GEOMETRIC") rd = 2. * (srad4[i] * srad4[k]);
            else if (radrule == "CUBIC-MEAN") {
                real rad4i2 = std::pow(rad4[i],2);
                real rad4k2 = std::pow(rad4[k],2);
                real rad4i3 = std::pow(rad4[i],3);
                real rad4k3 = std::pow(rad4[k],3);
                rd = 2. * (rad4i3+rad4k3) / (rad4i2+rad4k2);
            }
            else rd = rad4[i] + rad4[k];
            radmin4[kk][ii] = rd;
            radmin4[ii][kk] = rd;
        }
    }

    // use combination rules to set pairwise 1-4 well depths
    for (int ii = 0; ii < nlist; ii++) {
        i = list[ii];
        for (int kk = ii; kk < nlist; kk++) {
            k = list[kk];
            if (epsrule == "MMFF94") {
                ep = 0.;
                if (nn[i]!=0. and nn[k]!=0. and radmin4[kk][ii]!=0.) {
                    ep = 181.16*g[i]*g[k]*alph[i]*alph[k] / ((std::sqrt(alph[i]/nn[i])+std::sqrt(alph[k]/nn[k]))*std::pow(radmin4[kk][ii],6));
                }
                if (i == k) eps4[i] = ep;
            }
            else if (eps4[i]==0. and eps4[k]==0.) ep = 0.;
            else if (epsrule == "ARITHMETIC") ep = 0.5 * (eps4[i] + eps4[k]);
            else if (epsrule == "GEOMETRIC") ep = seps4[i] * seps4[k];
            else if (epsrule == "HARMONIC") ep = 2. * (eps4[i]*eps4[k]) / (eps4[i]+eps4[k]);
            else if (epsrule == "HHG") ep = 4. * (eps4[i]*eps4[k]) / std::pow((seps4[i]+seps4[k]),2);
            else if (epsrule == "W-H") ep = 2. * (seps4[i]*seps4[k]) * std::pow((rad4[i]*rad4[k]),3) / (std::pow(rad4[i],6)+std::pow(rad4[k],6));
            else ep = seps4[i] * seps4[k];
            epsilon4[kk][ii] = ep;
            epsilon4[ii][kk] = ep;
        }
    }

    // use reduced values for MMFF donor-acceptor pairs
    if (forcefield == "MMFF94") {
        for (int ii = 0; ii < nlist; ii++) {
            i = list[ii];
            for (int kk = ii; kk < nlist; kk++) {
                k = list[kk];
                if ((da[i]=='D' and da[k]=='A') or (da[i]=='A' and da[k]=='D')) {
                    epsilon[kk][ii] = epsilon[kk][ii] * 0.5;
                    epsilon[ii][kk] = epsilon[ii][kk] * 0.5;
                    radmin[kk][ii] = radmin[kk][ii] * 0.8;
                    radmin[ii][kk] = radmin[ii][kk] * 0.8;
                    epsilon4[kk][ii] = epsilon4[kk][ii] * 0.5;
                    epsilon4[ii][kk] = epsilon4[ii][kk] * 0.5;
                    radmin4[kk][ii] = radmin4[kk][ii] * 0.8;
                    radmin4[ii][kk] = radmin4[ii][kk] * 0.8;
                }
            }
        }
    }

    // vdw reduction factor information for each individual atom
    for (int i = 0; i < n; i++) {
        ired[i] = i;
        kred[i] = 0.;
        if (vdwindex == "TYPE") {
            kred[i] = reduct[type[i]];
        }
        else {
            kred[i] = reduct[atomClass[i]];
        }
        if (n12[i]==1 and kred[i]!=0.) {
            ired[i] = i12[i][0];
        }
    }

    // apply radii and well depths for special atom class pairs
    for (int i = 0; i < maxnvp; i++) {
        if (kvpr[i] == blank) break;
        ia = std::stoi(kvpr[i].substr(0,4));
        ib = std::stoi(kvpr[i].substr(4,4));
        ia--;
        ib--;
        if (rad[ia] == 0.) rad[ia] = 0.001;
        if (rad[ib] == 0.) rad[ib] = 0.001;
        ia = mvdw[ia];
        ib = mvdw[ib];
        if (ia!=-1 and ib!=-1) {
            if (radtyp == "SIGMA") radpr[i] = twosix * radpr[i];
            radmin[ib][ia] = radpr[i];
            radmin[ia][ib] = radpr[i];
            epsilon[ib][ia] = std::abs(epspr[i]);
            epsilon[ia][ib] = std::abs(epspr[i]);
            radmin4[ib][ia] = radpr[i];
            radmin4[ia][ib] = radpr[i];
            epsilon4[ib][ia] = std::abs(epspr[i]);
            epsilon4[ia][ib] = std::abs(epspr[i]);
        }
    }

    // set radii and well depths for hydrogen bonding pairs
    if (vdwtyp == "MM3-HBOND") {
        for (int i = 0; i < nlist; i++) {
            for (int k = 0; k < nlist; k++) {
                radhbnd[i][k] = 0.;
                epshbnd[i][k] = 0.;
            }
        }
        for (int i = 0; i < maxnhb; i++) {
            if (khb[i] == blank) break;
            ia = std::stoi(khb[i].substr(0,4));
            ib = std::stoi(khb[i].substr(4,4));
            ia--;
            ib--;
            if (rad[ia] == 0.) rad[ia] = 0.001;
            if (rad[ib] == 0.) rad[ib] = 0.001;
            ia = mvdw[ia];
            ib = mvdw[ib];
            if (ia!=-1 and ib!=-1) {
                if (radtyp == "SIGMA") radhb[i] = twosix * radhb[i];
                radhbnd[ib][ia] = radhb[i];
                radhbnd[ia][ib] = radhb[i];
                epshbnd[ib][ia] = std::abs(epshb[i]);
                epshbnd[ia][ib] = std::abs(epshb[i]);
            }
        }
    }

    // perform deallocation of some local arrays
    list.resize(0);
    srad.resize(0);
    srad4.resize(0);
    seps.resize(0);
    seps4.resize(0);

    // set coefficients for Gaussian fit to eps=1 and radmin=1
    if (vdwtyp == "GAUSSIAN") {
        real twosix2 = twosix * twosix;
        if (gausstyp == "LJ-4") {
            ngauss = 4;
            igauss[0][0] = 846706.7;
            igauss[0][1] = 15.464405 * twosix2;
            igauss[1][0] = 2713.651;
            igauss[1][1] = 7.346875 * twosix2;
            igauss[2][0] = -9.699172;
            igauss[2][1] = 1.8503725 * twosix2;
            igauss[3][0] = -0.715442;
            igauss[3][1] = 0.639621 * twosix2;
        }
        else if (gausstyp == "LJ-2") {
            ngauss = 2;
            igauss[0][0] = 14487.1;
            igauss[0][1] = 9.05148 * twosix2;
            igauss[1][0] = -5.55338;
            igauss[1][1] = 1.22536 * twosix2;
        }
        else if (gausstyp == "MM3-2") {
            ngauss = 2;
            igauss[0][0] = 2438.886;
            igauss[0][1] = 9.342616;
            igauss[1][0] = -6.197368;
            igauss[1][1] = 1.564486;
        }
        else if (gausstyp == "MM2-2") {
            ngauss = 2;
            igauss[0][0] = 3423.562;
            igauss[0][1] = 9.692821;
            igauss[1][0] = -6.50376;
            igauss[1][1] = 1.585344;
        }
        else if (gausstyp == "IN-PLACE") {
            ngauss = 2;
            igauss[0][0] = 500.;
            igauss[0][1] = 6.143;
            igauss[1][0] = -18.831;
            igauss[1][1] = 2.209;
        }
    }

    // remove zero-sized atoms from the list of vdw sites
    nvdw = 0;
    for (int i = 0; i < n; i++) {
        if (jvdw[i] != -1) {
            k = atomClass[i];
            if (vdwindex == "TYPE") k = type[i];
            if (rad[k] != 0.) {
                ivdw[nvdw] = i;
                nvdw++;
            }
        }
    }

    // turn off the van der Waals potential if it is not used
    if (nvdw == 0) use_vdw = false;
}
}
