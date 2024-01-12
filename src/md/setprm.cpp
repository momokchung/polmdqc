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
    if (bcon.size() != 0) bcon.resize(0);
    bcon.resize(maxnb);
    if (bcon5.size() != 0) bcon5.resize(0);
    bcon5.resize(maxnb5);
    if (bcon4.size() != 0) bcon4.resize(0);
    bcon4.resize(maxnb4);
    if (bcon3.size() != 0) bcon3.resize(0);
    bcon3.resize(maxnb3);
    if (blen.size() != 0) blen.resize(0);
    blen.resize(maxnb);
    if (blen5.size() != 0) blen5.resize(0);
    blen5.resize(maxnb5);
    if (blen4.size() != 0) blen4.resize(0);
    blen4.resize(maxnb4);
    if (blen3.size() != 0) blen3.resize(0);
    blen3.resize(maxnb3);
    if (dlen.size() != 0) dlen.resize(0);
    dlen.resize(maxnel);
    if (kb.size() != 0) kb.resize(0);
    kb.resize(maxnb);
    if (kb5.size() != 0) kb5.resize(0);
    kb5.resize(maxnb5);
    if (kb4.size() != 0) kb4.resize(0);
    kb4.resize(maxnb4);
    if (kb3.size() != 0) kb3.resize(0);
    kb3.resize(maxnb3);
    if (kel.size() != 0) kel.resize(0);
    kel.resize(maxnel);

    // allocate bond angle bend forcefield parameters
    if (acon.size() != 0) acon.resize(0);
    acon.resize(maxna);
    if (acon5.size() != 0) acon5.resize(0);
    acon5.resize(maxna5);
    if (acon4.size() != 0) acon4.resize(0);
    acon4.resize(maxna4);
    if (acon3.size() != 0) acon3.resize(0);
    acon3.resize(maxna3);
    if (aconp.size() != 0) aconp.resize(0);
    aconp.resize(maxnap);
    if (aconf.size() != 0) aconf.resize(0);
    aconf.resize(maxnaf);
    if (ang.size() != 0) ang.resize(0);
    ang.resize(maxna, std::vector<real>(3));
    if (ang5.size() != 0) ang5.resize(0);
    ang5.resize(maxna5, std::vector<real>(3));
    if (ang4.size() != 0) ang4.resize(0);
    ang4.resize(maxna4, std::vector<real>(3));
    if (ang3.size() != 0) ang3.resize(0);
    ang3.resize(maxna3, std::vector<real>(3));
    if (angp.size() != 0) angp.resize(0);
    angp.resize(maxnap, std::vector<real>(2));
    if (angf.size() != 0) angf.resize(0);
    angf.resize(maxnaf, std::vector<real>(2));
    if (ka.size() != 0) ka.resize(0);
    ka.resize(maxna);
    if (ka5.size() != 0) ka5.resize(0);
    ka5.resize(maxna5);
    if (ka4.size() != 0) ka4.resize(0);
    ka4.resize(maxna4);
    if (ka3.size() != 0) ka3.resize(0);
    ka3.resize(maxna3);
    if (kap.size() != 0) kap.resize(0);
    kap.resize(maxnap);
    if (kaf.size() != 0) kaf.resize(0);
    kaf.resize(maxnaf);

    // allocate stretch-bend forcefield parameters
    if (stbn.size() != 0) stbn.resize(0);
    stbn.resize(maxnsb, std::vector<real>(2));
    if (ksb.size() != 0) ksb.resize(0);
    ksb.resize(maxnsb);

    // allocate Urey-Bradley term forcefield parameters
    if (ucon.size() != 0) ucon.resize(0);
    ucon.resize(maxnu);
    if (dst13.size() != 0) dst13.resize(0);
    dst13.resize(maxnu);
    if (ku.size() != 0) ku.resize(0);
    ku.resize(maxnu);

    // allocate out-of-plane bend forcefield parameters
    if (opbn.size() != 0) opbn.resize(0);
    opbn.resize(maxnopb);
    if (kopb.size() != 0) kopb.resize(0);
    kopb.resize(maxnopb);

    // allocate out-of-plane distance forcefield parameters
    if (opds.size() != 0) opds.resize(0);
    opds.resize(maxnopd);
    if (kopd.size() != 0) kopd.resize(0);
    kopd.resize(maxnopd);

    // allocate improper dihedral forcefield parameters
    if (dcon.size() != 0) dcon.resize(0);
    dcon.resize(maxndi);
    if (tdi.size() != 0) tdi.resize(0);
    tdi.resize(maxndi);
    if (kdi.size() != 0) kdi.resize(0);
    kdi.resize(maxndi);

    // allocate improper torsion forcefield parameters
    if (ti1.size() != 0) ti1.resize(0);
    ti1.resize(maxnti, std::vector<real>(2));
    if (ti2.size() != 0) ti2.resize(0);
    ti2.resize(maxnti, std::vector<real>(2));
    if (ti3.size() != 0) ti3.resize(0);
    ti3.resize(maxnti, std::vector<real>(2));
    if (kti.size() != 0) kti.resize(0);
    kti.resize(maxnti);

    // allocate torsion angle forcefield parameters
    if (t1.size() != 0) t1.resize(0);
    t1.resize(maxnt, std::vector<real>(2));
    if (t2.size() != 0) t2.resize(0);
    t2.resize(maxnt, std::vector<real>(2));
    if (t3.size() != 0) t3.resize(0);
    t3.resize(maxnt, std::vector<real>(2));
    if (t4.size() != 0) t4.resize(0);
    t4.resize(maxnt, std::vector<real>(2));
    if (t5.size() != 0) t5.resize(0);
    t5.resize(maxnt, std::vector<real>(2));
    if (t6.size() != 0) t6.resize(0);
    t6.resize(maxnt, std::vector<real>(2));
    if (t15.size() != 0) t15.resize(0);
    t15.resize(maxnt5, std::vector<real>(2));
    if (t25.size() != 0) t25.resize(0);
    t25.resize(maxnt5, std::vector<real>(2));
    if (t35.size() != 0) t35.resize(0);
    t35.resize(maxnt5, std::vector<real>(2));
    if (t45.size() != 0) t45.resize(0);
    t45.resize(maxnt5, std::vector<real>(2));
    if (t55.size() != 0) t55.resize(0);
    t55.resize(maxnt5, std::vector<real>(2));
    if (t65.size() != 0) t65.resize(0);
    t65.resize(maxnt5, std::vector<real>(2));
    if (t14.size() != 0) t14.resize(0);
    t14.resize(maxnt4, std::vector<real>(2));
    if (t24.size() != 0) t24.resize(0);
    t24.resize(maxnt4, std::vector<real>(2));
    if (t34.size() != 0) t34.resize(0);
    t34.resize(maxnt4, std::vector<real>(2));
    if (t44.size() != 0) t44.resize(0);
    t44.resize(maxnt4, std::vector<real>(2));
    if (t54.size() != 0) t54.resize(0);
    t54.resize(maxnt4, std::vector<real>(2));
    if (t64.size() != 0) t64.resize(0);
    t64.resize(maxnt4, std::vector<real>(2));
    if (kt.size() != 0) kt.resize(0);
    kt.resize(maxnt);
    if (kt5.size() != 0) kt5.resize(0);
    kt5.resize(maxnt5);
    if (kt4.size() != 0) kt4.resize(0);
    kt4.resize(maxnt4);

    // allocate pi-system torsion forcefield parameters
    if (ptcon.size() != 0) ptcon.resize(0);
    ptcon.resize(maxnpt);
    if (kpt.size() != 0) kpt.resize(0);
    kpt.resize(maxnpt);

    // allocate stretch-torsion forcefield parameters
    if (btcon.size() != 0) btcon.resize(0);
    btcon.resize(maxnbt, std::vector<real>(9));
    if (kbt.size() != 0) kbt.resize(0);
    kbt.resize(maxnbt);

    // allocate angle-torsion forcefield parameters
    if (atcon.size() != 0) atcon.resize(0);
    atcon.resize(maxnat, std::vector<real>(6));
    if (kat.size() != 0) kat.resize(0);
    kat.resize(maxnat);

    // allocate torsion-torsion forcefield parameters
    if (tnx.size() != 0) tnx.resize(0);
    tnx.resize(maxntt);
    if (tny.size() != 0) tny.resize(0);
    tny.resize(maxntt);
    if (ttx.size() != 0) ttx.resize(0);
    ttx.resize(maxntt, std::vector<real>(maxtgrd));
    if (tty.size() != 0) tty.resize(0);
    tty.resize(maxntt, std::vector<real>(maxtgrd));
    if (tbf.size() != 0) tbf.resize(0);
    tbf.resize(maxntt, std::vector<real>(maxtgrd2));
    if (tbx.size() != 0) tbx.resize(0);
    tbx.resize(maxntt, std::vector<real>(maxtgrd2));
    if (tby.size() != 0) tby.resize(0);
    tby.resize(maxntt, std::vector<real>(maxtgrd2));
    if (tbxy.size() != 0) tbxy.resize(0);
    tbxy.resize(maxntt, std::vector<real>(maxtgrd2));
    if (ktt.size() != 0) ktt.resize(0);
    ktt.resize(maxntt);

    // allocate special vdw term forcefield parameters
    if (radpr.size() != 0) radpr.resize(0);
    radpr.resize(maxnvp);
    if (epspr.size() != 0) epspr.resize(0);
    epspr.resize(maxnvp);
    if (kvpr.size() != 0) kvpr.resize(0);
    kvpr.resize(maxnvp);

    // allocate H-bonding term forcefield parameters
    if (radhb.size() != 0) radhb.resize(0);
    radhb.resize(maxnhb);
    if (epshb.size() != 0) epshb.resize(0);
    epshb.resize(maxnhb);
    if (khb.size() != 0) khb.resize(0);
    khb.resize(maxnhb);

    // allocate bond dipole forcefield parameters
    if (dpl.size() != 0) dpl.resize(0);
    dpl.resize(maxnd);
    if (dpl5.size() != 0) dpl5.resize(0);
    dpl5.resize(maxnd5);
    if (dpl4.size() != 0) dpl4.resize(0);
    dpl4.resize(maxnd4);
    if (dpl3.size() != 0) dpl3.resize(0);
    dpl3.resize(maxnd3);
    if (pos.size() != 0) pos.resize(0);
    pos.resize(maxnd);
    if (pos5.size() != 0) pos5.resize(0);
    pos5.resize(maxnd5);
    if (pos4.size() != 0) pos4.resize(0);
    pos4.resize(maxnd4);
    if (pos3.size() != 0) pos3.resize(0);
    pos3.resize(maxnd3);
    if (kd.size() != 0) kd.resize(0);
    kd.resize(maxnd);
    if (kd5.size() != 0) kd5.resize(0);
    kd5.resize(maxnd5);
    if (kd4.size() != 0) kd4.resize(0);
    kd4.resize(maxnd4);
    if (kd3.size() != 0) kd3.resize(0);
    kd3.resize(maxnd3);

    // allocate atomic multipole forcefield parameters
    if (multip.size() != 0) multip.resize(0);
    multip.resize(maxnmp, std::vector<real>(13));
    if (mpaxis.size() != 0) mpaxis.resize(0);
    mpaxis.resize(maxnmp);
    if (kmp.size() != 0) kmp.resize(0);
    kmp.resize(maxnmp);

    // allocate special Thole forcefield parameters
    if (thlpr.size() != 0) thlpr.resize(0);
    thlpr.resize(maxnpp);
    if (thdpr.size() != 0) thdpr.resize(0);
    thdpr.resize(maxnpp);
    if (kppr.size() != 0) kppr.resize(0);
    kppr.resize(maxnpp);

    // allocate charge flux term forcefield parameters
    if (cflb.size() != 0) cflb.resize(0);
    cflb.resize(maxncfb);
    if (cfla.size() != 0) cfla.resize(0);
    cfla.resize(maxncfa, std::vector<real>(2));
    if (cflab.size() != 0) cflab.resize(0);
    cflab.resize(maxncfa, std::vector<real>(2));
    if (kcfb.size() != 0) kcfb.resize(0);
    kcfb.resize(maxncfb);
    if (kcfa.size() != 0) kcfa.resize(0);
    kcfa.resize(maxncfa);

    // allocate pisystem orbital forcefield parameters
    if (sslope.size() != 0) sslope.resize(0);
    sslope.resize(maxnpi);
    if (sslope5.size() != 0) sslope5.resize(0);
    sslope5.resize(maxnpi5);
    if (sslope4.size() != 0) sslope4.resize(0);
    sslope4.resize(maxnpi4);
    if (tslope.size() != 0) tslope.resize(0);
    tslope.resize(maxnpi);
    if (tslope5.size() != 0) tslope5.resize(0);
    tslope5.resize(maxnpi5);
    if (tslope4.size() != 0) tslope4.resize(0);
    tslope4.resize(maxnpi4);
    if (kpi.size() != 0) kpi.resize(0);
    kpi.resize(maxnpi);
    if (kpi5.size() != 0) kpi5.resize(0);
    kpi5.resize(maxnpi5);
    if (kpi4.size() != 0) kpi4.resize(0);
    kpi4.resize(maxnpi4);
}
}
