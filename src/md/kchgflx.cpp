// Author: Moses KJ Chung
// Year:   2024

#include "sizes.h"
#include "angbnd.h"
#include "atmlst.h"
#include "atomid.h"
#include "atoms.h"
#include "bndstr.h"
#include "cflux.h"
#include "couple.h"
#include "gettext.h"
#include "inform.h"
#include "kangs.h"
#include "kbonds.h"
#include "kcflux.h"
#include "kchgflx.h"
#include "keys.h"
#include "numeral.h"
#include "potent.h"
#include "upcase.h"
#include "usage.h"
#include <algorithm>
#include <sstream>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  kchgflx  --  charge flux parameter assignment  //
//                                                 //
/////////////////////////////////////////////////////

// "kchgflx" assigns bond stretch and angle bend charge flux
// correction values and processes any new or changed values
// for these parameters

void kchgflx()
{
    int ia,ib,ic;
    int ita,itb,itc;
    int na,nb;
    int size,next;
    double cfb;
    double cfa1,cfa2;
    double cfb1,cfb2;
    bool header;
    std::string pa,pb,pc;
    std::string blank,pt;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // process keywords containing bond charge flux parameters
    blank = "";
    size = 4;
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "BNDCFLUX") {
            ia = 0;
            ib = 0;
            cfb = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> cfb;
            if (std::min(ia,ib) <= 0) continue;
            if (header and !silent) {
                header = false;
                printf("\n Additional Bond Charge Flux Parameters :");
                printf("\n\n     Atom Classes                   K(CFB)\n\n");
            }
            if (!silent) {
                printf("      %4d%4d             %15.6f\n", ia,ib,cfb);
            }
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
                pt = pa + pb;
            }
            else {
                pt = pb + pa;
            }
            for (int j = 0; j < maxncfb; j++) {
                if (kcfb[j]==blank or kcfb[j]==pt) {
                    kcfb[j] = pt;
                    if (ia < ib) {
                        cflb[j] = cfb;
                    }
                    else if (ib < ia) {
                        cflb[j] = -cfb;
                    }
                    else {
                        cflb[j] = 0.;
                        printf("\n KCHGFLX  --  Bond Charge Flux for Identical Classes Set to Zero\n");
                    }
                    break;
                }
            }
        }
    }

    // process keywords containing angle charge flux parameters
    blank = "";
    size = 4;
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "ANGCFLUX") {
            ia = 0;
            ib = 0;
            ic = 0;
            cfa1 = 0.;
            cfa2 = 0.;
            cfb1 = 0.;
            cfb2 = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> ic >> cfa1 >> cfa2 >> cfb1 >> cfb2;
            if (std::min({ia,ib,ic}) <= 0) continue;
            if (header and !silent) {
                header = false;
                printf("\n Additional Angle Charge Flux Parameters :");
                printf("\n\n     Atom Classes          K(CFA1)       K(CFA2)       K(CFB1)       K(CFB2)\n\n");
            }
            if (!silent) {
                printf("    %4d%4d%4d    %14.6f%14.6f%14.6f%14.6f\n", ia,ib,ic,cfa1,cfa2,cfb1,cfb2);
            }
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            pc = numeral(ic,size);
            if (ia <= ic) {
                pt = pa + pb + pc;
            }
            else {
                pt = pc + pb + pa;
            }
            for (int j = 0; j < maxncfa; j++) {
                if (kcfa[j]==blank or kcfa[j]==pt) {
                    kcfa[j] = pt;
                    if (ia <= ic) {
                        cfla[j][0] = cfa1;
                        cfla[j][1] = cfa2;
                        cflab[j][0] = cfb1;
                        cflab[j][1] = cfb2;
                    }
                    else {
                        cfla[j][0] = cfa2;
                        cfla[j][1] = cfa1;
                        cflab[j][0] = cfb2;
                        cflab[j][1] = cfb1;
                    }
                    break;
                }
            }
        }
    }

    // determine the total number of forcefield parameters
    nb = maxncfb;
    for (int i = maxncfb-1; i >= 0; i--) {
        if (kcfb[i] == blank) nb = i;
    }
    na = maxncfa;
    for (int i = maxncfa-1; i >= 0; i--) {
        if (kcfa[i] == blank) na = i;
    }

    // perform dynamic allocation of some global arrays
    bflx.allocate(nbond);
    aflx.allocate(nangle);
    abflx.allocate(nangle);

    // assign bond charge flux parameters for each bond
    nbflx = 0;
    for (int i = 0; i < nbond; i++) {
        ia = ibnd[i][0];
        ib = ibnd[i][1];
        ita = atomClass[ia];
        itb = atomClass[ib];
        size = 4;
        pa = numeral(ita+1,size);
        pb = numeral(itb+1,size);
        if (ita <= itb) {
            pt = pa + pb;
        }
        else {
            pt = pb + pa;
        }
        bflx[i] = 0.;
        for (int j = 0; j < nb; j++) {
            if (kcfb[j] == pt) {
                nbflx += 1;
                if (ita <= itb) {
                    bflx[i] = cflb[j];
                }
                else {
                    bflx[i] = -cflb[j];
                }
            }
        }
    }

    // assign angle charge flux parameters for each angle
    naflx = 0;
    for (int i = 0; i < nangle; i++) {
        ia = iang[i][0];
        ib = iang[i][1];
        ic = iang[i][2];
        ita = atomClass[ia];
        itb = atomClass[ib];
        itc = atomClass[ic];
        pa = numeral(ita+1,size);
        pb = numeral(itb+1,size);
        pc = numeral(itc+1,size);
        if (ita <= itc) {
            pt = pa + pb + pc;
        }
        else {
            pt = pc + pb + pa;
        }
        aflx[i][0] = 0.;
        aflx[i][1] = 0.;
        abflx[i][0] = 0.;
        abflx[i][1] = 0.;
        for (int j = 0; j < na; j++) {
            if (kcfa[j] == pt) {
                naflx += 1;
                if (ita <= itc) {
                    aflx[i][0] = cfla[j][0];
                    aflx[i][1] = cfla[j][1];
                    abflx[i][0] = cflab[j][0];
                    abflx[i][1] = cflab[j][1];
                }
                else {
                    aflx[i][0] = cfla[j][1];
                    aflx[i][1] = cfla[j][0];
                    abflx[i][0] = cflab[j][1];
                    abflx[i][1] = cflab[j][0];
                }
            }
        }
    }

    // process keywords with bond charge flux parameters for specific bonds
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "BNDCFLUX") {
            ia = 0;
            ib = 0;
            cfb = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> cfb;
            if (std::min(ia,ib) >= 0) continue;
            ia = std::abs(ia);
            ib = std::abs(ib);
            if (header and !silent) {
                header = false;
                printf("\n Additional Bond Charge Flux Parameters for Specific Bonds :");
                printf("\n\n        Atoms                       K(CFB)\n\n");
            }
            if (!silent) {
                printf("      %4d%4d             %15.6f\n", ia,ib,cfb);
            }
            for (int j = 0; j < nbond; j++) {
                ita = ibnd[j][0] + 1;
                itb = ibnd[j][1] + 1;
                if ((ia==ita and ib==itb) or (ib==itb and ib==ita)) {
                    if (ia < ib) {
                        bflx[j] = cfb;
                    }
                    else if (ib < ia) {
                        bflx[j] = -cfb;
                    }
                    else {
                        bflx[j] = 0.;
                        printf("\n KCHGFLX  --  Bond Charge Flux for Identical Classes Set to Zero\n");
                    }
                    break;
                }
            }
        }
    }

    // process keywords with angle charge flux parameters for specific angles
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "ANGCFLUX") {
            ia = 0;
            ib = 0;
            ic = 0;
            cfa1 = 0.;
            cfa2 = 0.;
            cfb1 = 0.;
            cfb2 = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> ic >> cfa1 >> cfa2 >> cfb1 >> cfb2;
            if (std::min({ia,ib,ic}) >= 0) continue;
            ia = std::abs(ia);
            ib = std::abs(ib);
            ic = std::abs(ic);
            if (header and !silent) {
                header = false;
                printf("\n Additional Angle Charge Flux Parameters for Specific Angles :");
                printf("\n\n        Atoms              K(CFA1)       K(CFA2)       K(CFB1)       K(CFB2)\n\n");
            }
            if (!silent) {
                printf("    %4d%4d%4d    %14.6f%14.6f%14.6f%14.6f\n", ia,ib,ic,cfa1,cfa2,cfb1,cfb2);
            }
            for (int j = 0; j < nangle; j++) {
                ita = iang[j][0] + 1;
                itb = iang[j][1] + 1;
                itc = iang[j][2] + 1;
                if (ib == itb) {
                    if ((ia==ita and ic==itc) or (ia==itc and ic==ita)) {
                        if (ia <= ic) {
                            aflx[i][0] = cfa1;
                            aflx[i][1] = cfa2;
                            abflx[i][0] = cfb1;
                            abflx[i][1] = cfb2;
                        }
                        else {
                            aflx[i][0] = cfa2;
                            aflx[i][1] = cfa1;
                            abflx[i][0] = cfb2;
                            abflx[i][1] = cfb1;
                        }
                        break;
                    }
                }
            }
        }
    }

    // turn off bond and angle charge flux term if not used
    if (nbflx==0 and naflx==0) use_chgflx = false;
    if (!use_charge and !use_mpole and !use_polar) use_chgflx = false;
}
}
