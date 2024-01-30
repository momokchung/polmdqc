// Author: Moses KJ Chung
// Year:   2024

#include "angbnd.h"
#include "angpot.h"
#include "atomid.h"
#include "atoms.h"
#include "chkring.h"
#include "couple.h"
#include "fields.h"
#include "gettext.h"
#include "inform.h"
#include "kangle.h"
#include "kangs.h"
#include "keys.h"
#include "numeral.h"
#include "potent.h"
#include "upcase.h"
#include "usage.h"
#include <algorithm>
#include <sstream>

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  kangle  --  angle bend parameter assignment  //
//                                               //
///////////////////////////////////////////////////

// "kangle" assigns the force constants and ideal angles for
// the bond angles; also processes new or changed parameters

void kangle()
{
    int ia,ib,ic;
    int ita,itb,itc;
    int na,nap,naf;
    int na3,na4,na5;
    int jen,ih,nh;
    int next,size;
    int minat,iring;
    real fc,an,pr;
    real an1,an2,an3;
    bool header,done;
    bool use_ring;
    std::string pa,pb,pc;
    std::string label;
    std::string blank,pt;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // process keywords containing angle bending parameters
    blank = "";
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        iring = -1;
        if (keyword == "ANGLE") iring = 0;
        if (keyword == "ANGLE5") iring = 5;
        if (keyword == "ANGLE4") iring = 4;
        if (keyword == "ANGLE3") iring = 3;
        if (iring >= 0) {
            ia = 0;
            ib = 0;
            ic = 0;
            fc = 0.;
            an1 = 0.;
            an2 = 0.;
            an3 = 0.;
            jen = 0;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> ic >> fc >> an1 >> an2 >> an3;
            if (std::min({ia,ib,ic}) <= 0) continue;
            if (an2!=0. or an3!=0.) jen = 1;
            if (!silent) {
                if (header) {
                    header = false;
                    printf("\n Additional Angle Bending Parameters :");
                    printf("\n\n     Atom Classes             K(B)          Angle\n\n");
                }
                if (iring == 0) {
                    if (jen == 0) {
                        printf("    %4d%4d%4d   %15.3f%15.3f\n", ia,ib,ic,fc,an1);
                    }
                    else if (an1 != 0.) {
                        printf("    %4d%4d%4d   %15.3f%15.3f   0-H's\n", ia,ib,ic,fc,an1);
                    }
                    if (an2 != 0.) {
                        printf("    %4d%4d%4d   %15.3f%15.3f   1-H's\n", ia,ib,ic,fc,an2);
                    }
                    if (an3 != 0.) {
                        printf("    %4d%4d%4d   %15.3f%15.3f   2-H's\n", ia,ib,ic,fc,an3);
                    }
                }
                else {
                    if (iring == 5) label = "5-Ring";
                    if (iring == 4) label = "4-Ring";
                    if (iring == 3) label = "3-Ring";
                    if (jen == 0) {
                        printf("    %4d%4d%4d   %15.3f%15.3f   %6s\n", ia,ib,ic,fc,an1,label.c_str());
                    }
                    else if (an1 != 0.) {
                        printf("    %4d%4d%4d   %15.3f%15.3f   %6s   0-H's\n", ia,ib,ic,fc,an1,label.c_str());
                    }
                    if (an2 != 0.) {
                        printf("    %4d%4d%4d   %15.3f%15.3f   %6s   1-H's\n", ia,ib,ic,fc,an2,label.c_str());
                    }
                    if (an3 != 0.) {
                        printf("    %4d%4d%4d   %15.3f%15.3f   %6s   2-H's\n", ia,ib,ic,fc,an3,label.c_str());
                    }
                }
            }
            size = 4;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            pc = numeral(ic,size);
            if (ia <= ic) {
                pt = pa + pb + pc;
            }
            else {
                pt = pc + pb + pa;
            }
            if (iring == 0) {
                for (int j = 0; j < maxna; j++) {
                    if (ka[j]==blank or ka[j]==pt) {
                        ka[j] = pt;
                        acon[j] = fc;
                        ang[j][0] = an1;
                        ang[j][1] = an2;
                        ang[j][2] = an3;
                        goto label1_end;
                    }
                }
                printf("\n KANGLE  --  Too many Bond Angle Bending Parameters\n");
                informAbort = true;
            }
            else if (iring == 5) {
                for (int j = 0; j < maxna5; j++) {
                    if (ka5[j]==blank or ka5[j]==pt) {
                        ka5[j] = pt;
                        acon5[j] = fc;
                        ang5[j][0] = an1;
                        ang5[j][1] = an2;
                        ang5[j][2] = an3;
                        goto label1_end;
                    }
                }
                printf("\n KANGLE  --  Too many 5-Ring Angle Bending Parameters\n");
                informAbort = true;
            }
            else if (iring == 4) {
                for (int j = 0; j < maxna4; j++) {
                    if (ka4[j]==blank or ka4[j]==pt) {
                        ka4[j] = pt;
                        acon4[j] = fc;
                        ang4[j][0] = an1;
                        ang4[j][1] = an2;
                        ang4[j][2] = an3;
                        goto label1_end;
                    }
                }
                printf("\n KANGLE  --  Too many 4-Ring Angle Bending Parameters\n");
                informAbort = true;
            }
            else if (iring == 3) {
                for (int j = 0; j < maxna3; j++) {
                    if (ka3[j]==blank or ka3[j]==pt) {
                        ka3[j] = pt;
                        acon3[j] = fc;
                        ang3[j][0] = an1;
                        ang3[j][1] = an2;
                        ang3[j][2] = an3;
                        goto label1_end;
                    }
                }
                printf("\n KANGLE  --  Too many 3-Ring Angle Bending Parameters\n");
                informAbort = true;
            }
            label1_end:;
        }
    }

    // process keywords containing in-plane angle bending parameters
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "ANGLEP") {
            ia = 0;
            ib = 0;
            ic = 0;
            fc = 0.;
            an1 = 0.;
            an2 = 0.;
            jen = 0;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> ic >> fc >> an1 >> an2;
            if (std::min({ia,ib,ic}) <= 0) continue;
            if (an2 != 0.) jen = 1;
            if (!silent) {
                if (header) {
                    header = false;
                    printf("\n Additional In-Plane Angle Bending Parameters :");
                    printf("\n\n     Atom Classes             K(B)          Angle\n\n");
                }
                if (jen == 0) {
                    printf("    %4d%4d%4d   %15.3f%15.3f\n", ia,ib,ic,fc,an1);
                }
                else if (an1 != 0.) {
                    printf("    %4d%4d%4d   %15.3f%15.3f   0-H's\n", ia,ib,ic,fc,an1);
                }
                if (an2 != 0.) {
                    printf("    %4d%4d%4d   %15.3f%15.3f   1-H's\n", ia,ib,ic,fc,an2);
                }
            }
            size = 4;
            pa = numeral (ia,size);
            pb = numeral (ib,size);
            pc = numeral (ic,size);
            if (ia <= ic) {
                pt = pa + pb + pc;
            }
            else {
                pt = pc + pb + pa;
            }
            for (int j = 0; j < maxnap; j++) {
                if (kap[j]==blank or kap[j]==pt) {
                    kap[j] = pt;
                    aconp[j] = fc;
                    angp[j][0] = an1;
                    angp[j][1] = an2;
                    goto label_260;;
                }
            }
            printf("\n KANGLE  --  Too many In-Plane Angle Bending Parameters\n");
            informAbort = true;
            label_260:;
        }
    }

    // process keywords containing Fourier angle bending parameters
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "ANGLEF") {
            ia = 0;
            ib = 0;
            ic = 0;
            fc = 0.;
            an = 0.;
            pr = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> ic >> fc >> an >> pr;
            if (std::min({ia,ib,ic}) <= 0) continue;
            if (!silent) {
                if (header) {
                    header = false;
                    printf("\n Additional Fourier Angle Bending Parameters :");
                    printf("\n\n     Atom Classes             K(B)          Shift         Period\n\n");
                }
                printf("    %4d%4d%4d   %15.3f%15.3f%15.3f\n", ia,ib,ic,fc,an,pr);
            }
            size = 4;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            pc = numeral(ic,size);
            if (ia <= ic) {
                pt = pa + pb + pc;
            }
            else {
                pt = pc + pb + pa;
            }
            for (int j = 0; j < maxnaf; j++) {
                if (kaf[j]==blank or kaf[j]==pt) {
                    kaf[j] = pt;
                    aconf[j] = fc;
                    angf[j][0] = an;
                    angf[j][1] = pr;
                    goto label_310;
                }
            }
            printf("\n KANGLE  --  Too many Fourier Angle Bending Parameters\n");
            informAbort = true;
            label_310:;
        }
    }

    // determine the total number of forcefield parameters
    na = maxna;
    na5 = maxna5;
    na4 = maxna4;
    na3 = maxna3;
    nap = maxnap;
    naf = maxnaf;
    for (int i = maxna-1; i >= 0; i--) {
        if (ka[i] == blank) na = i;
    }
    for (int i = maxna5-1; i >= 0; i--) {
        if (ka5[i] == blank) na5 = i;
    }
    for (int i = maxna4-1; i >= 0; i--) {
        if (ka4[i] == blank) na4 = i;
    }
    for (int i = maxna3-1; i >= 0; i--) {
        if (ka3[i] == blank) na3 = i;
    }
    for (int i = maxnap-1; i >= 0; i--) {
        if (kap[i] == blank) nap = i;
    }
    for (int i = maxnaf-1; i >= 0; i--) {
        if (kaf[i] == blank) naf = i;
    }
    use_ring = false;
    if (std::min({na5,na4,na3}) != 0) use_ring = true;

    // set generic parameters for use with any number of hydrogens
    for (int i = 0; i < na; i++) {
        if (ang[i][1]==0. and ang[i][2]==0.) {
            ang[i][1] = ang[i][0];
            ang[i][2] = ang[i][0];
        }
    }
    for (int i = 0; i < na5; i++) {
        if (ang5[i][1]==0. and ang5[i][2]==0.) {
            ang5[i][1] = ang5[i][0];
            ang5[i][2] = ang5[i][0];
        }
    }
    for (int i = 0; i < na4; i++) {
        if (ang4[i][1]==0. and ang4[i][2]==0.) {
            ang4[i][1] = ang4[i][0];
            ang4[i][2] = ang4[i][0];
        }
    }
    for (int i = 0; i < na3; i++) {
        if (ang3[i][1]==0. and ang3[i][2]==0.) {
            ang3[i][1] = ang3[i][0];
            ang3[i][2] = ang3[i][0];
        }
    }
    for (int i = 0; i < nap; i++) {
        if (angp[i][1] == 0.) {
            angp[i][1] = angp[i][0];
        }
    }

    // perform dynamic allocation of some global arrays
    ak.allocate(nangle);
    anat.allocate(nangle);
    afld.allocate(nangle);
    angtyp.allocate(nangle);

    // assign ideal bond angle and force constant for each angle
    header = true;
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
        if (ita <= itc) {
            pt = pa + pb + pc;
        }
        else {
            pt = pc + pb + pa;
        }
        ak[i] = 0.;
        anat[i] = 0.;
        afld[i] = 0.;
        angtyp[i] = "HARMONIC";
        done = false;

        // count number of non-angle hydrogens on the central atom
        nh = 1;
        for (int j = 0; j < n12[ib] ; j++) {
            ih = i12[ib][j];
            if (ih!=ia and ih!=ic and atomic[ih]==1) nh += 1;
        }

        // make a check for bond angles contained inside small rings
        iring = 0;
        if (use_ring) {
            chkring(iring,ia,ib,ic,-1);
            if (iring == 6) iring = 0;
            if (iring==5 and na5==0) iring = 0;
            if (iring==4 and na4==0) iring = 0;
            if (iring==3 and na3==0) iring = 0;
        }

        // assign angle bending parameters for bond angles
        if (iring == 0) {
            for (int j = 0; j < na; j++) {
                if (ka[j]==pt and ang[j][nh]!=0.) {
                    ak[i] = acon[j];
                    anat[i] = ang[j][nh];
                    done = true;
                    goto label_320;
                }
            }
        }

        // assign bending parameters for 5-membered ring angles
        else if (iring == 5) {
            for (int j = 0; j < na5; j++) {
                if (ka5[j]==pt and ang5[j][nh]!=0.) {
                    ak[i] = acon5[j];
                    anat[i] = ang5[j][nh];
                    done = true;
                    goto label_320;
                }
            }
        }

        // assign bending parameters for 4-membered ring angles
        else if (iring == 4) {
            for (int j = 0; j < na4; j++) {
                if (ka4[j]==pt and ang4[j][nh]!=0.) {
                    ak[i] = acon4[j];
                    anat[i] = ang4[j][nh];
                    done = true;
                    goto label_320;
                }
            }
        }

        // assign bending parameters for 3-membered ring angles
        else if (iring == 3) {
            for (int j = 0; j < na3; j++) {
                if (ka3[j]==pt and ang3[j][nh]!=0.) {
                    ak[i] = acon3[j];
                    anat[i] = ang3[j][nh];
                    done = true;
                    goto label_320;
                }
            }
        }

        // assign in-plane angle bending parameters for bond angles
        if (!done and n12[ib]==3) {
            for (int j = 0; j < nap; j++) {
                if (kap[j]==pt and angp[j][nh]!=0.) {
                    ak[i] = aconp[j];
                    anat[i] = angp[j][nh];
                    angtyp[i] = "IN-PLANE";
                    done = true;
                    goto label_320;
                }
            }
        }

        // assign Fourier angle bending parameters for bond angles
        if (!done) {
            for (int j = 0; j < naf; j++) {
                if (kaf[j] == pt) {
                    ak[i] = aconf[j];
                    anat[i] = angf[j][0];
                    afld[i] = angf[j][1];
                    angtyp[i] = "FOURIER";
                    done = true;
                    goto label_320;
                }
            }
        }

        // warning if suitable angle bending parameter not found
        label_320:
        minat = std::min({atomic[ia],atomic[ib],atomic[ic]});
        if (minat == 0) done = true;
        if (use_angle and !done) {
            if (use[ia+1] or use[ib+1] or use[ic+1]) informAbort = true;
            if (header) {
                header = false;
                printf("\n Undefined Angle Bending Parameters :");
                printf("\n\n Type                  Atom Names                   Atom Classes\n\n");
            }
            label = "Angle ";
            if (iring == 5) label = "5-Ring";
            if (iring == 4) label = "4-Ring";
            if (iring == 3) label = "3-Ring";
            printf(" %6s     %6d-%-3s%6d-%-3s%6d-%-3s       %5d%5d%5d\n",
                label.c_str(),ia+1,name[ia].c_str(),ib+1,name[ib].c_str(),ic+1,name[ic].c_str(),ita+1,itb+1,itc+1);
        }
    }

    // process keywords containing angle specific parameters
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "ANGLE") {
            ia = 0;
            ib = 0;
            ic = 0;
            fc = 0.;
            an1 = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> ic >> fc >> an1;
            if (std::min({ia,ib,ic}) < 0) {
                ia = std::abs(ia);
                ib = std::abs(ib);
                ic = std::abs(ic);
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Angle Parameters for Specific Angles :");
                    printf("\n\n        Atoms                 K(B)          Angle\n\n");
                }
                if (!silent) {
                    printf("    %4d%4d%4d   %15.3f%15.3f\n", ia,ib,ic,fc,an1);
                }
                for (int j = 0; j < nangle; j++) {
                    ita = iang[j][0];
                    itb = iang[j][1];
                    itc = iang[j][2];
                    ia--;
                    ib--;
                    ic--;
                    if (ib == itb) {
                        if ((ia==ita and ic==itc) or (ia==itc and ic==ita)) {
                            ak[j] = fc;
                            anat[j] = an1;
                            angtyp[j] = "HARMONIC";
                            goto label_380;
                        }
                    }
                }
            }
            label_380:;
        }
    }

    // process keywords containing in-plane angle specific parameters
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "ANGLEP") {
            ia = 0;
            ib = 0;
            ic = 0;
            fc = 0.;
            an1 = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> ic >> fc >> an1;
            if (std::min({ia,ib,ic}) < 0) {
                ia = std::abs(ia);
                ib = std::abs(ib);
                ic = std::abs(ic);
                if (header and !silent) {
                    header = false;
                    printf("\n Additional In-Plane Angle Parameters for Specific Angles :");
                    printf("\n\n        Atoms                 K(B)          Angle\n\n");
                }
                if (!silent) {
                    printf("    %4d%4d%4d   %15.3f%15.3f\n", ia,ib,ic,fc,an1);
                }
                if (ia > ic) {
                    ita = ia;
                    ia = ic;
                    ic = ita;
                }
                ia--;
                ib--;
                ic--;
                for (int j = 0; j < nangle; j++) {
                    if (ia==iang[j][0] and ib==iang[j][1] and ic==iang[j][2]) {
                        ak[j] = fc;
                        anat[j] = an1;
                        angtyp[j] = "IN-PLANE";
                        goto label_420;
                    }
                }
            }
            label_420:;
        }
    }

    // process keywords containing Fourier angle specific parameters
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "ANGLEF") {
            ia = 0;
            ib = 0;
            ic = 0;
            fc = 0.;
            an = 0.;
            pr = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia >> ib >> ic >> fc >> an >> pr;
            if (std::min({ia,ib,ic}) < 0) {
                ia = std::abs(ia);
                ib = std::abs(ib);
                ic = std::abs(ic);
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Fourier Angle Parameters for Specific Angles :");
                    printf("\n\n        Atoms                 K(B)          Shift         Period\n\n");
                }
                if (!silent) {
                    printf("    %4d%4d%4d   %15.3f%15.3f%15.3f\n", ia,ib,ic,fc,an,pr);
                }
                if (ia > ic) {
                    ita = ia;
                    ia = ic;
                    ic = ita;
                }
                ia--;
                ib--;
                ic--;
                for (int j = 0; j < nangle; j++) {
                    if (ia==iang[j][0] and ib==iang[j][1] and ic==iang[j][2]) {
                        ak[j] = fc;
                        anat[j] = an;
                        afld[j] = pr;
                        angtyp[j] = "FOURIER";
                        goto label_460;
                    }
                }
            }
            label_460:;
        }
    }

    // turn off the angle bending potential if it is not used
    if (nangle == 0) use_angle = false;
}
}
