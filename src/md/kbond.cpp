// Author: Moses KJ Chung
// Year:   2024

#include "angbnd.h"
#include "atmlst.h"
#include "atomid.h"
#include "atoms.h"
#include "bndstr.h"
#include "chkring.h"
#include "couple.h"
#include "fields.h"
#include "gettext.h"
#include "inform.h"
#include "kbond.h"
#include "kbonds.h"
#include "keys.h"
#include "numeral.h"
#include "potent.h"
#include "tors.h"
#include "upcase.h"
#include "usage.h"
#include <algorithm>
#include <sstream>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  kbond  --  bond stretch parameter assignment  //
//                                                //
////////////////////////////////////////////////////

// "kbond" assigns a force constant and ideal bond length
// to each bond in the structure and processes any new or
// changed parameter values

void kbond()
{
    int ia,ib,ita,itb;
    int nb,nb5,nb4,nb3;
    int size,next;
    int minat,iring;
    real fc,bd;
    bool header,done;
    bool use_ring;
    std::string pa,pb;
    std::string label;
    std::string blank,pt;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // process keywords containing bond stretch parameters
    blank = "";
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        iring = -1;
        if (keyword == "BOND") iring = 0;
        if (keyword == "BOND5") iring = 5;
        if (keyword == "BOND4") iring = 4;
        if (keyword == "BOND3") iring = 3;
        if (iring >= 0) {
            ia = 0;
            ib = 0;
            fc = 0.;
            bd = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> fc >> bd;
            if (std::min(ia,ib) <= 0) goto label1_end;
            if (!silent) {
                if (header) {
                    header = false;
                    printf("\n Additional Bond Stretching Parameters :");
                    printf("\n\n     Atom Classes             K(S)         Length\n\n");
                }
                if (iring == 0) {
                    printf("      %4d%4d     %15.3f%15.4f\n",ia,ib,fc,bd);
                }
                else {
                    if (iring == 5) label = "5-Ring";
                    if (iring == 4) label = "4-Ring";
                    if (iring == 3) label = "3-Ring";
                    printf("      %4d%4d     %15.3f%15.4f   %s\n",ia,ib,fc,bd,label.c_str());
                }
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
            if (iring == 0) {
                for (int j = 0; j < maxnb; j++) {
                    if (kb[j]==blank or kb[j]==pt) {
                        kb[j] = pt;
                        bcon[j] = fc;
                        blen[j] = bd;
                        goto label1_end;
                    }
                }
                printf("\n KBOND  --  Too many Bond Stretching Parameters\n");
                informAbort = true;
            }
            else if (iring == 5) {
                for (int j = 0; j < maxnb5; j++) {
                    if (kb5[j]==blank or kb5[j]==pt) {
                        kb5[j] = pt;
                        bcon5[j] = fc;
                        blen5[j] = bd;
                        goto label1_end;
                    }
                }
                printf("\n KBOND  --  Too many 5-Ring Stretching Parameters\n");
                informAbort = true;
            }
            else if (iring == 4) {
                for (int j = 0; j < maxnb4; j++) {
                    if (kb4[j]==blank or kb4[j]==pt) {
                        kb4[j] = pt;
                        bcon4[j] = fc;
                        blen4[j] = bd;
                        goto label1_end;
                    }
                }
                printf("\n KBOND  --  Too many 4-Ring Stretching Parameters\n");
                informAbort = true;
            }
            else if (iring == 3) {
                for (int j = 0; j < maxnb3; j++) {
                    if (kb3[j]==blank or kb3[j]==pt) {
                        kb3[j] = pt;
                        bcon3[j] = fc;
                        blen3[j] = bd;
                        goto label1_end;
                    }
                }
                printf("\n KBOND  --  Too many 3-Ring Stretching Parameters\n");
                informAbort = true;
            }
            label1_end:
            continue;
        }
    }

    // determine the total number of forcefield parameters
    nb = maxnb;
    nb5 = maxnb5;
    nb4 = maxnb4;
    nb3 = maxnb3;
    for (int i = maxnb-1; i >= 0; i--) {
        if (kb[i] == blank) nb = i;
    }
    for (int i = maxnb5-1; i >= 0; i--) {
        if (kb5[i] == blank) nb5 = i;
    }
    for (int i = maxnb4-1; i >= 0; i--) {
        if (kb4[i] == blank) nb4 = i;
    }
    for (int i = maxnb3-1; i >= 0; i--) {
        if (kb3[i] == blank) nb3 = i;
    }
    use_ring = false;
    if (std::min({nb5,nb4,nb3}) != 0) use_ring = true;

    // perform dynamic allocation of some global arrays
    bk.allocate(nbond);
    bl.allocate(nbond);

    // assign ideal bond length and force constant for each bond
    header = true;
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
        bk[i] = 0.;
        bl[i] = 0.;
        done = false;

        // make a check for bonds contained inside small rings
        iring = 0;
        if (use_ring) {
            chkring(iring,ia,ib,-1,-1);
            if (iring == 6) iring = 0;
            if (iring==5 and nb5==0) iring = 0;
            if (iring==4 and nb4==0) iring = 0;
            if (iring==3 and nb3==0) iring = 0;
        }

        // assign bond stretching parameters for each bond
        if (iring == 0) {
            for (int j = 0; j < nb; j++) {
                if (kb[j] == pt) {
                    bk[i] = bcon[j];
                    bl[i] = blen[j];
                    done = true;
                    goto label2_end;
                }
            }
        }

        // assign stretching parameters for 5-membered ring bonds
        else if (iring == 5) {
            for (int j = 0; j < nb5; j++) {
                if (kb5[j] == pt) {
                    bk[i] = bcon5[j];
                    bl[i] = blen5[j];
                    done = true;
                    goto label2_end;
                }
            }
        }

        // assign stretching parameters for 4-membered ring bonds
        else if (iring == 4) {
            for (int j = 0; j < nb4; j++) {
                if (kb4[j] == pt) {
                    bk[i] = bcon4[j];
                    bl[i] = blen4[j];
                    done = true;
                    goto label2_end;
                }
            }
        }

        // assign stretching parameters for 3-membered ring bonds
        else if (iring == 3) {
            for (int j = 0; j < nb3; j++) {
                if (kb3[j] == pt) {
                    bk[i] = bcon3[j];
                    bl[i] = blen3[j];
                    done = true;
                    goto label2_end;
                }
            }
        }

        // warning if suitable bond stretching parameter not found
        label2_end:;
        minat = std::min(atomic[ia],atomic[ib]);
        if (minat == 0) done = true;
        if (use_bond and !done) {
            if (use[ia+1] or use[ib+1]) informAbort = true;
            if (header) {
                header = false;
                printf("\n Undefined Bond Stretching Parameters :");
                printf("\n\n Type             Atom Names           Atom Classes\n\n");
            }
            label = "Bond  ";
            if (iring == 5) label = "5-Ring";
            if (iring == 4) label = "4-Ring";
            if (iring == 3) label = "3-Ring";
            printf(" %6s     %6d-%-3s%6d-%-3s       %5d%5d\n", label.c_str(),ia+1,name[ia].c_str(),ib+1,name[ib].c_str(),ita+1,itb+1);
        }
    }

    // process keywords containing bond specific parameters
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "BOND") {
            ia = 0;
            ib = 0;
            fc = 0.;
            bd = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> fc >> bd;
            if (std::min(ia,ib) < 0) {
                ia = std::abs(ia);
                ib = std::abs(ib);
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Bond Parameters for Specific Bonds :");
                    printf("\n\n        Atoms                 K(S)         Length\n\n");
                }
                if (!silent) {
                    printf("      %4d%4d     %15.3f%15.4f\n",ia,ib,fc,bd);
                }
                for (int j = 0; j < nbond; j++) {
                    ita = ibnd[j][0] + 1;
                    itb = ibnd[j][1] + 1;
                    if ((ia==ita and ib==itb) or (ib==itb and ib==ita)) {
                        bk[j] = fc;
                        bl[j] = bd;
                        goto label_200;
                    }
                }
            }
            label_200:;
        }
    }

    // check for electronegativity bond length corrections
    keneg();

    // turn off the bond stretch potential if it is not used
    if (nbond == 0) use_bond = false;
}

//////////////////////////////////////////////////////
//                                                  //
//  keneg  --  assign electronegativity parameters  //
//                                                  //
//////////////////////////////////////////////////////

// "keneg" applies primary and secondary electronegativity bond
// length corrections to applicable bond parameters

// note this version does not scale multiple corrections to the
// same bond by increasing powers of 0.62 as in MM3

void keneg()
{
    int ia,ib,ic,id,m,nel;
    int ita,itb,itc,itd;
    int size,next;
    real dl,factor;
    bool header;
    std::string pa,pb,pc,pd;
    std::string blank;
    std::string pt,pt1,pt2;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // process keywords containing electronegativity parameters
    blank = "";
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "ELECTNEG") {
            ia = 0;
            ib = 0;
            ic = 0;
            dl = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> ic >> dl;
            if (!silent) {
                if (header) {
                    header = false;
                    printf("\n Additional Electronegativity Parameters :");
                    printf("\n\n     Atom Classes                  dLength\n\n");
                }
                printf("    %4d%4d%4d              %12.4f\n",ia,ib,ic,dl);
            }
            size = 4;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            pc = numeral(ic,size);
            pt = pa + pb + pc;
            for (int j = 0; j < maxnel; j++) {
                if (kel[j]==blank or kel[j]==pt) {
                    kel[j] = pt;
                    dlen[j] = dl;
                    goto label_50;
                }
            }
            printf("\n KENEG  --  Too many Electronegativity Parameters\n");
            informAbort = true;
            label_50:;
        }
    }

    // determine the total number of forcefield parameters
    nel = maxnel;
    for (int i = 0; i < maxnel; i++) {
        if (kel[i] == blank) {
            nel = i;
            break;
        }
    }

    // check angles for primary electronegativity corrections
    if (nel != 0) {
        for (int i = 0; i < nangle; i++) {
            ia = iang[i][0];
            ib = iang[i][1];
            ic = iang[i][2];
            ita = atomClass[ia];
            itb = atomClass[ib];
            itc = atomClass[ic];
            size = 4;
            pa = numeral(ita+1,size);
            pb = numeral(itb+1,size);
            pc = numeral(itc+1,size);
            pt1 = pa + pb + pc;
            pt2 = pc + pb + pa;

            // search the parameter set for a match to either bond
            for (int j = 0; j < nel; j++) {
                if (kel[j] == pt1) {
                    for (int k = 0; k < n12[ia]; k++) {
                        if (i12[ia][k] == ib) {
                            m = bndlist[ia][k];
                            bl[m] += dlen[j];
                        }
                    }
                    goto label_70;
                }
                else if (kel[j] == pt2) {
                    for (int k = 0; k < n12[ic]; k++) {
                        if (i12[ic][k] == ib) {
                            m = bndlist[ic][k];
                            bl[m] += dlen[j];
                        }
                    }
                    goto label_70;
                }
            }
            label_70:;
        }

        // check torsions for secondary electronegativity corrections
        factor = 0.4;
        for (int i = 0; i < ntors; i++) {
            ia = itors[i][0];
            ib = itors[i][1];
            ic = itors[i][2];
            id = itors[i][3];
            ita = atomClass[ia];
            itb = atomClass[ib];
            itc = atomClass[ic];
            itd = atomClass[id];
            size = 4;
            pa = numeral(ita+1,size);
            pb = numeral(itb+1,size);
            pc = numeral(itc+1,size);
            pd = numeral(itd+1,size);
            pt1 = pa + pb + pd;
            pt2 = pd + pc + pa;

            // turn off electronegativity effect for attached hydrogens
            if (atomic[id] <= 1) pt1 = blank;
            if (atomic[ia] <= 1) pt2 = blank;

            // search the parameter set for a match to either bond
            for (int j = 0; j < nel; j++) {
                if (kel[j] == pt1) {
                    for (int k = 0; k < n12[ia]; k++) {
                        if (i12[ia][k] == ib) {
                            m = bndlist[ia][k];
                            bl[m] += factor*dlen[j];
                        }
                    }
                    goto label_80;
                }
                else if (kel[j] == pt2) {
                    for (int k = 0; k < n12[id]; k++) {
                        if (i12[id][k] == ic) {
                            m = bndlist[id][k];
                            bl[m] += factor*dlen[j];
                        }
                    }
                    goto label_80;
                }
            }
            label_80:;
        }
    }
}
}
