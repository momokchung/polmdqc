// Author: Moses KJ Chung
// Year:   2023

#include "atoms.h"
#include "couple.h"
#include "kangs.h"
#include "kantor.h"
#include "kbonds.h"
#include "kcflux.h"
#include "kdipol.h"
#include "keys.h"
#include "khbond.h"
#include "kiprop.h"
#include "kitors.h"
#include "kmulti.h"
#include "kopbnd.h"
#include "kopdst.h"
#include "korbs.h"
#include "kpitor.h"
#include "kpolpr.h"
#include "kstbnd.h"
#include "ksttor.h"
#include "ktorsn.h"
#include "ktrtor.h"
#include "kurybr.h"
#include "kvdwpr.h"
#include "params.h"
#include "restrn.h"
#include "setprm.h"
#include "upcase.h"
#include <string>
#include <sstream>

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  setprm  --  allocate force field parameters  //
//                                               //
///////////////////////////////////////////////////

// "setprm" counts and allocates memory space for force field
// parameter values involving multiple atom types or classes as
// found in the parameter file and keyfile

void setprm()
{
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // zero out the count of each force field parameter type
    maxnb = 0;
    maxnb5 = 0;
    maxnb4 = 0;
    maxnb3 = 0;
    maxnel = 0;
    maxna = 0;
    maxna5 = 0;
    maxna4 = 0;
    maxna3 = 0;
    maxnap = 0;
    maxnaf = 0;
    maxnsb = 0;
    maxnu = 0;
    maxnopb = 0;
    maxnopd = 0;
    maxndi = 0;
    maxnti = 0;
    maxnt = 0;
    maxnt5 = 0;
    maxnt4 = 0;
    maxnpt = 0;
    maxnbt = 0;
    maxnat = 0;
    maxntt = 0;
    maxnvp = 0;
    maxnhb = 0;
    maxnd = 0;
    maxnd5 = 0;
    maxnd4 = 0;
    maxnd3 = 0;
    maxnmp = 0;
    maxnpp = 0;
    maxncfb = 0;
    maxncfa = 0;
    maxnpi = 0;
    maxnpi5 = 0;
    maxnpi4 = 0;
    maxfix = 0;

    // find any parameter values found in the parameter file
    for (int i = 0; i < nprm; i++) {
        record = prmline[i];
        iss.clear();
        keyword = "";
        iss.str(record);
        iss >> keyword;
        upcase(keyword);
        if (keyword == "BOND") maxnb++;
        if (keyword == "MMFFBOND") maxnb++;
        if (keyword == "BOND5") maxnb5++;
        if (keyword == "BOND4") maxnb4++;
        if (keyword == "BOND3") maxnb3++;
        if (keyword == "ELECTNEG") maxnel++;
        if (keyword == "ANGLE") maxna++;
        if (keyword == "MMFFANGLE") maxna++;
        if (keyword == "ANGLE5") maxna5++;
        if (keyword == "ANGLE4") maxna4++;
        if (keyword == "ANGLE3") maxna3++;
        if (keyword == "ANGLEP") maxnap++;
        if (keyword == "ANGLEF") maxnaf++;
        if (keyword == "STRBND") maxnsb++;
        if (keyword == "UREYBRAD") maxnu++;
        if (keyword == "OPBEND") maxnopb++;
        if (keyword == "OPDIST") maxnopd++;
        if (keyword == "IMPROPER") maxndi++;
        if (keyword == "IMPTORS") maxnti++;
        if (keyword == "TORSION") maxnt++;
        if (keyword == "MMFFTORSION") maxnt++;
        if (keyword == "TORSION5") maxnt5++;
        if (keyword == "TORSION4") maxnt4++;
        if (keyword == "PITORS") maxnpt++;
        if (keyword == "STRTORS") maxnbt++;
        if (keyword == "ANGTORS") maxnat++;
        if (keyword == "TORTORS") maxntt++;
        if (keyword == "VDWPAIR") maxnvp++;
        if (keyword == "VDWPR") maxnvp++;
        if (keyword == "HBOND") maxnhb++;
        if (keyword == "DIPOLE") maxnd++;
        if (keyword == "DIPOLE5") maxnd5++;
        if (keyword == "DIPOLE4") maxnd4++;
        if (keyword == "DIPOLE3") maxnd3++;
        if (keyword == "MULTIPOLE") maxnmp++;
        if (keyword == "POLPAIR") maxnpp++;
        if (keyword == "BNDCFLUX") maxncfb++;
        if (keyword == "ANGCFLUX") maxncfa++;
        if (keyword == "PIBOND") maxnpi++;
        if (keyword == "PIBOND5") maxnpi5++;
        if (keyword == "PIBOND4") maxnpi4++;
        if (keyword.substr(0,9) == "RESTRAIN-") {
            maxfix++;
            if (keyword.substr(9) == "POSITION") {
                string = keyword.substr(18);
                int ia = 0;
                iss >> ia;
                if (ia>=-n and ia<=-1) {
                    ia = std::abs(ia);
                    int ib = 0;
                    iss >> ib;
                    ib = std::min(std::abs(ib), n);
                    maxfix += std::max(0,ib-ia);
                }
            }
        }
        if (keyword == "ENFORCE-CHIRALITY") {
            for (int k = 0; k < n; k++) {
                if (n12[k] == 4) maxfix++;
            }
        }
    }

    // find additional parameter values found in the keyfile
    for (int i = 0; i < nkey; i++) {
        record = keyline[i];
        iss.clear();
        keyword = "";
        iss.str(record);
        iss >> keyword;
        upcase(keyword);
        if (keyword == "BOND") maxnb++;
        if (keyword == "MMFFBOND") maxnb++;
        if (keyword == "BOND5") maxnb5++;
        if (keyword == "BOND4") maxnb4++;
        if (keyword == "BOND3") maxnb3++;
        if (keyword == "ELECTNEG") maxnel++;
        if (keyword == "ANGLE") maxna++;
        if (keyword == "MMFFANGLE") maxna++;
        if (keyword == "ANGLE5") maxna5++;
        if (keyword == "ANGLE4") maxna4++;
        if (keyword == "ANGLE3") maxna3++;
        if (keyword == "ANGLEP") maxnap++;
        if (keyword == "ANGLEF") maxnaf++;
        if (keyword == "STRBND") maxnsb++;
        if (keyword == "UREYBRAD") maxnu++;
        if (keyword == "OPBEND") maxnopb++;
        if (keyword == "OPDIST") maxnopd++;
        if (keyword == "IMPROPER") maxndi++;
        if (keyword == "IMPTORS") maxnti++;
        if (keyword == "TORSION") maxnt++;
        if (keyword == "MMFFTORSION") maxnt++;
        if (keyword == "TORSION5") maxnt5++;
        if (keyword == "TORSION4") maxnt4++;
        if (keyword == "PITORS") maxnpt++;
        if (keyword == "STRTORS") maxnbt++;
        if (keyword == "ANGTORS") maxnat++;
        if (keyword == "TORTORS") maxntt++;
        if (keyword == "VDWPAIR") maxnvp++;
        if (keyword == "VDWPR") maxnvp++;
        if (keyword == "HBOND") maxnhb++;
        if (keyword == "DIPOLE") maxnd++;
        if (keyword == "DIPOLE5") maxnd5++;
        if (keyword == "DIPOLE4") maxnd4++;
        if (keyword == "DIPOLE3") maxnd3++;
        if (keyword == "MULTIPOLE") maxnmp++;
        if (keyword == "POLPAIR") maxnpp++;
        if (keyword == "BNDCFLUX") maxncfb++;
        if (keyword == "ANGCFLUX") maxncfa++;
        if (keyword == "PIBOND") maxnpi++;
        if (keyword == "PIBOND5") maxnpi5++;
        if (keyword == "PIBOND4") maxnpi4++;
        if (keyword.substr(0,9) == "RESTRAIN-") {
            maxfix++;
            if (keyword.substr(9) == "POSITION") {
                string = keyword.substr(18);
                int ia = 0;
                iss >> ia;
                if (ia>=-n and ia<=-1) {
                    ia = std::abs(ia);
                    int ib = 0;
                    iss >> ib;
                    ib = std::min(std::abs(ib), n);
                    maxfix += std::max(0,ib-ia);
                }
            }
        }
        if (keyword == "ENFORCE-CHIRALITY") {
            for (int k = 0; k < n; k++) {
                if (n12[k] == 4) maxfix++;
            }
        }
    }

    // set the allocated memory for each parameter type
    maxnb = std::max(500,maxnb+100);
    maxnb5 = std::max(500,maxnb5+100);
    maxnb4 = std::max(500,maxnb4+100);
    maxnb3 = std::max(500,maxnb3+100);
    maxnel = std::max(500,maxnel+100);
    maxna = std::max(500,maxna+100);
    maxna5 = std::max(500,maxna5+100);
    maxna4 = std::max(500,maxna4+100);
    maxna3 = std::max(500,maxna3+100);
    maxnap = std::max(500,maxnap+100);
    maxnaf = std::max(500,maxnaf+100);
    maxnsb = std::max(500,maxnsb+100);
    maxnu = std::max(500,maxnu+100);
    maxnopb = std::max(500,maxnopb+100);
    maxnopd = std::max(500,maxnopd+100);
    maxndi = std::max(500,maxndi+100);
    maxnti = std::max(500,maxnti+100);
    maxnt = std::max(500,maxnt+100);
    maxnt5 = std::max(500,maxnt5+100);
    maxnt4 = std::max(500,maxnt4+100);
    maxnpt = std::max(500,maxnpt+100);
    maxnbt = std::max(500,maxnbt+100);
    maxnat = std::max(500,maxnat+100);
    maxntt = std::max(50,maxntt+10);
    maxnvp = std::max(500,maxnvp+100);
    maxnhb = std::max(500,maxnhb+100);
    maxnd = std::max(500,maxnd+100);
    maxnd5 = std::max(500,maxnd5+100);
    maxnd4 = std::max(500,maxnd4+100);
    maxnd3 = std::max(500,maxnd3+100);
    maxnmp = std::max(500,maxnmp+100);
    maxnpp = std::max(500,maxnpp+100);
    maxncfb = std::max(500,maxncfb+100);
    maxncfa = std::max(500,maxncfa+100);
    maxnpi = std::max(500,maxnpi+100);
    maxnpi5 = std::max(500,maxnpi5+100);
    maxnpi4 = std::max(500,maxnpi4+100);
    maxfix = std::max(500,maxfix+100);

    // allocate bond stretching forcefield parameters
    bcon.allocate(maxnb);
    bcon5.allocate(maxnb5);
    bcon4.allocate(maxnb4);
    bcon3.allocate(maxnb3);
    blen.allocate(maxnb);
    blen5.allocate(maxnb5);
    blen4.allocate(maxnb4);
    blen3.allocate(maxnb3);
    dlen.allocate(maxnel);
    kb.allocate(maxnb);
    kb5.allocate(maxnb5);
    kb4.allocate(maxnb4);
    kb3.allocate(maxnb3);
    kel.allocate(maxnel);

    // allocate bond angle bend forcefield parameters
    acon.allocate(maxna);
    acon5.allocate(maxna5);
    acon4.allocate(maxna4);
    acon3.allocate(maxna3);
    aconp.allocate(maxnap);
    aconf.allocate(maxnaf);
    ang.allocate(maxna);
    ang5.allocate(maxna5);
    ang4.allocate(maxna4);
    ang3.allocate(maxna3);
    angp.allocate(maxnap);
    angf.allocate(maxnaf);
    ka.allocate(maxna);
    ka5.allocate(maxna5);
    ka4.allocate(maxna4);
    ka3.allocate(maxna3);
    kap.allocate(maxnap);
    kaf.allocate(maxnaf);

    // allocate stretch-bend forcefield parameters
    stbn.allocate(maxnsb);
    ksb.allocate(maxnsb);

    // allocate Urey-Bradley term forcefield parameters
    ucon.allocate(maxnu);
    dst13.allocate(maxnu);
    ku.allocate(maxnu);

    // allocate out-of-plane bend forcefield parameters
    opbn.allocate(maxnopb);
    kopb.allocate(maxnopb);

    // allocate out-of-plane distance forcefield parameters
    opds.allocate(maxnopd);
    kopd.allocate(maxnopd);

    // allocate improper dihedral forcefield parameters
    dcon.allocate(maxndi);
    tdi.allocate(maxndi);
    kdi.allocate(maxndi);

    // allocate improper torsion forcefield parameters
    ti1.allocate(maxnti);
    ti2.allocate(maxnti);
    ti3.allocate(maxnti);
    kti.allocate(maxnti);

    // allocate torsion angle forcefield parameters
    t1.allocate(maxnt);
    t2.allocate(maxnt);
    t3.allocate(maxnt);
    t4.allocate(maxnt);
    t5.allocate(maxnt);
    t6.allocate(maxnt);
    t15.allocate(maxnt5);
    t25.allocate(maxnt5);
    t35.allocate(maxnt5);
    t45.allocate(maxnt5);
    t55.allocate(maxnt5);
    t65.allocate(maxnt5);
    t14.allocate(maxnt4);
    t24.allocate(maxnt4);
    t34.allocate(maxnt4);
    t44.allocate(maxnt4);
    t54.allocate(maxnt4);
    t64.allocate(maxnt4);
    kt.allocate(maxnt);
    kt5.allocate(maxnt5);
    kt4.allocate(maxnt4);

    // allocate pi-system torsion forcefield parameters
    ptcon.allocate(maxnpt);
    kpt.allocate(maxnpt);

    // allocate stretch-torsion forcefield parameters
    btcon.allocate(maxnbt);
    kbt.allocate(maxnbt);

    // allocate angle-torsion forcefield parameters
    atcon.allocate(maxnat);
    kat.allocate(maxnat);

    // allocate torsion-torsion forcefield parameters
    tnx.allocate(maxntt);
    tny.allocate(maxntt);
    ttx.allocate(maxntt);
    tty.allocate(maxntt);
    tbf.allocate(maxntt);
    tbx.allocate(maxntt);
    tby.allocate(maxntt);
    tbxy.allocate(maxntt);
    ktt.allocate(maxntt);

    // allocate special vdw term forcefield parameters
    radpr.allocate(maxnvp);
    epspr.allocate(maxnvp);
    kvpr.allocate(maxnvp);

    // allocate H-bonding term forcefield parameters
    radhb.allocate(maxnhb);
    epshb.allocate(maxnhb);
    khb.allocate(maxnhb);

    // allocate bond dipole forcefield parameters
    dpl.allocate(maxnd);
    dpl5.allocate(maxnd5);
    dpl4.allocate(maxnd4);
    dpl3.allocate(maxnd3);
    pos.allocate(maxnd);
    pos5.allocate(maxnd5);
    pos4.allocate(maxnd4);
    pos3.allocate(maxnd3);
    kd.allocate(maxnd);
    kd5.allocate(maxnd5);
    kd4.allocate(maxnd4);
    kd3.allocate(maxnd3);

    // allocate atomic multipole forcefield parameters
    multip.allocate(maxnmp);
    mpaxis.allocate(maxnmp);
    kmp.allocate(maxnmp);

    // allocate special Thole forcefield parameters
    thlpr.allocate(maxnpp);
    thdpr.allocate(maxnpp);
    kppr.allocate(maxnpp);

    // allocate charge flux term forcefield parameters
    cflb.allocate(maxncfb);
    cfla.allocate(maxncfa);
    cflab.allocate(maxncfa);
    kcfb.allocate(maxncfb);
    kcfa.allocate(maxncfa);

    // allocate pisystem orbital forcefield parameters
    sslope.allocate(maxnpi);
    sslope5.allocate(maxnpi5);
    sslope4.allocate(maxnpi4);
    tslope.allocate(maxnpi);
    tslope5.allocate(maxnpi5);
    tslope4.allocate(maxnpi4);
    kpi.allocate(maxnpi);
    kpi5.allocate(maxnpi5);
    kpi4.allocate(maxnpi4);
}
}
