///////////////////////////////////////////////////////
//                                                   //
//  setprm.cpp  --  allocate force field parameters  //
//                                                   //
///////////////////////////////////////////////////////

// "setprm" counts and allocates memory space for force field
// parameter values involving multiple atom types or classes as
// found in the parameter file and keyfile


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
    bcon.resize(maxnb);
    bcon5.resize(maxnb5);
    bcon4.resize(maxnb4);
    bcon3.resize(maxnb3);
    blen.resize(maxnb);
    blen5.resize(maxnb5);
    blen4.resize(maxnb4);
    blen3.resize(maxnb3);
    dlen.resize(maxnel);
    kb.resize(maxnb);
    kb5.resize(maxnb5);
    kb4.resize(maxnb4);
    kb3.resize(maxnb3);
    kel.resize(maxnel);

    // allocate bond angle bend forcefield parameters
    acon.resize(maxna);
    acon5.resize(maxna5);
    acon4.resize(maxna4);
    acon3.resize(maxna3);
    aconp.resize(maxnap);
    aconf.resize(maxnaf);
    ang.resize(maxna, std::vector<double>(3));
    ang5.resize(maxna5, std::vector<double>(3));
    ang4.resize(maxna4, std::vector<double>(3));
    ang3.resize(maxna3, std::vector<double>(3));
    angp.resize(maxnap, std::vector<double>(2));
    angf.resize(maxnaf, std::vector<double>(2));
    ka.resize(maxna);
    ka5.resize(maxna5);
    ka4.resize(maxna4);
    ka3.resize(maxna3);
    kap.resize(maxnap);
    kaf.resize(maxnaf);

    // allocate stretch-bend forcefield parameters
    stbn.resize(maxnsb, std::vector<double>(2));
    ksb.resize(maxnsb);

    // allocate Urey-Bradley term forcefield parameters
    ucon.resize(maxnu);
    dst13.resize(maxnu);
    ku.resize(maxnu);

    // allocate out-of-plane bend forcefield parameters
    opbn.resize(maxnopb);
    kopb.resize(maxnopb);

    // allocate out-of-plane distance forcefield parameters
    opds.resize(maxnopd);
    kopd.resize(maxnopd);

    // allocate improper dihedral forcefield parameters
    dcon.resize(maxndi);
    tdi.resize(maxndi);
    kdi.resize(maxndi);

    // allocate improper torsion forcefield parameters
    ti1.resize(maxnti, std::vector<double>(2));
    ti2.resize(maxnti, std::vector<double>(2));
    ti3.resize(maxnti, std::vector<double>(2));
    kti.resize(maxnti);

    // allocate torsion angle forcefield parameters
    t1.resize(maxnt, std::vector<double>(2));
    t2.resize(maxnt, std::vector<double>(2));
    t3.resize(maxnt, std::vector<double>(2));
    t4.resize(maxnt, std::vector<double>(2));
    t5.resize(maxnt, std::vector<double>(2));
    t6.resize(maxnt, std::vector<double>(2));
    t15.resize(maxnt5, std::vector<double>(2));
    t25.resize(maxnt5, std::vector<double>(2));
    t35.resize(maxnt5, std::vector<double>(2));
    t45.resize(maxnt5, std::vector<double>(2));
    t55.resize(maxnt5, std::vector<double>(2));
    t65.resize(maxnt5, std::vector<double>(2));
    t14.resize(maxnt4, std::vector<double>(2));
    t24.resize(maxnt4, std::vector<double>(2));
    t34.resize(maxnt4, std::vector<double>(2));
    t44.resize(maxnt4, std::vector<double>(2));
    t54.resize(maxnt4, std::vector<double>(2));
    t64.resize(maxnt4, std::vector<double>(2));
    kt.resize(maxnt);
    kt5.resize(maxnt5);
    kt4.resize(maxnt4);

    // allocate pi-system torsion forcefield parameters
    ptcon.resize(maxnpt);
    kpt.resize(maxnpt);

    // allocate stretch-torsion forcefield parameters
    btcon.resize(maxnbt, std::vector<double>(9));
    kbt.resize(maxnbt);

    // allocate angle-torsion forcefield parameters
    atcon.resize(maxnat, std::vector<double>(6));
    kat.resize(maxnat);

    // allocate torsion-torsion forcefield parameters
    tnx.resize(maxntt);
    tny.resize(maxntt);
    ttx.resize(maxntt, std::vector<double>(maxtgrd));
    tty.resize(maxntt, std::vector<double>(maxtgrd));
    tbf.resize(maxntt, std::vector<double>(maxtgrd2));
    tbx.resize(maxntt, std::vector<double>(maxtgrd2));
    tby.resize(maxntt, std::vector<double>(maxtgrd2));
    tbxy.resize(maxntt, std::vector<double>(maxtgrd2));
    ktt.resize(maxntt);

    // allocate special vdw term forcefield parameters
    radpr.resize(maxnvp);
    epspr.resize(maxnvp);
    kvpr.resize(maxnvp);

    // allocate H-bonding term forcefield parameters
    radhb.resize(maxnhb);
    epshb.resize(maxnhb);
    khb.resize(maxnhb);

    // allocate bond dipole forcefield parameters
    dpl.resize(maxnd);
    dpl5.resize(maxnd5);
    dpl4.resize(maxnd4);
    dpl3.resize(maxnd3);
    pos.resize(maxnd);
    pos5.resize(maxnd5);
    pos4.resize(maxnd4);
    pos3.resize(maxnd3);
    kd.resize(maxnd);
    kd5.resize(maxnd5);
    kd4.resize(maxnd4);
    kd3.resize(maxnd3);

    // allocate atomic multipole forcefield parameters
    multip.resize(maxnmp, std::vector<double>(13));
    mpaxis.resize(maxnmp);
    kmp.resize(maxnmp);

    // allocate special Thole forcefield parameters
    thlpr.resize(maxnpp);
    thdpr.resize(maxnpp);
    kppr.resize(maxnpp);

    // allocate charge flux term forcefield parameters
    cflb.resize(maxncfb);
    cfla.resize(maxncfa, std::vector<double>(2));
    cflab.resize(maxncfa, std::vector<double>(2));
    kcfb.resize(maxncfb);
    kcfa.resize(maxncfa);

    // allocate pisystem orbital forcefield parameters
    sslope.resize(maxnpi);
    sslope5.resize(maxnpi5);
    sslope4.resize(maxnpi4);
    tslope.resize(maxnpi);
    tslope5.resize(maxnpi5);
    tslope4.resize(maxnpi4);
    kpi.resize(maxnpi);
    kpi5.resize(maxnpi5);
    kpi4.resize(maxnpi4);
}
