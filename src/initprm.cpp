//////////////////////////////////////////////////////////
//                                                      //
//  initprm.cpp  --  initialize force field parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// "initprm" completely initializes a force field by setting all
// parameters to zero and using defaults for control values


#include "angpot.h"
#include "bndpot.h"
#include "chgpot.h"
#include "ctrpot.h"
#include "dsppot.h"
#include "expol.h"
#include "extfld.h"
#include "fields.h"
#include "ielscf.h"
#include "initprm.h"
#include "kanang.h"
#include "kangs.h"
#include "kantor.h"
#include "katoms.h"
#include "kbonds.h"
#include "kcflux.h"
#include "kchrge.h"
#include "kcpen.h"
#include "kctrn.h"
#include "kdipol.h"
#include "kdsp.h"
#include "kexpl.h"
#include "khbond.h"
#include "kiprop.h"
#include "kitors.h"
#include "kmulti.h"
#include "kopbnd.h"
#include "kopdst.h"
#include "korbs.h"
#include "kpitor.h"
#include "kpolpr.h"
#include "kpolr.h"
#include "krepl.h"
#include "ksolut.h"
#include "kstbnd.h"
#include "ksttor.h"
#include "ktorsn.h"
#include "ktrtor.h"
#include "kurybr.h"
#include "kvdws.h"
#include "kvdwpr.h"
#include "mathConst.h"
#include "mplpot.h"
#include "polpot.h"
#include "reppot.h"
#include "rxnpot.h"
#include "solpot.h"
#include "sizes.h"
#include "solute.h"
#include "urypot.h"
#include "torpot.h"
#include "units.h"
#include "uprior.h"
#include "vdwpot.h"
#include <cmath>
#include <string>

void initmmff();

void initprm()
{
    // initialize strings of parameter atom types and classes
    for (int i = 0; i < maxnvp; i++) {
        kvpr[i] = "";
    }
    for (int i = 0; i < maxnhb; i++) {
        khb[i] = "";
    }
    for (int i = 0; i < maxnb; i++) {
        kb[i] = "";
    }
    for (int i = 0; i < maxnb5; i++) {
        kb5[i] = "";
    }
    for (int i = 0; i < maxnb4; i++) {
        kb4[i] = "";
    }
    for (int i = 0; i < maxnb3; i++) {
        kb3[i] = "";
    }
    for (int i = 0; i < maxnel; i++) {
        kel[i] = "";
    }
    for (int i = 0; i < maxna; i++) {
        ka[i] = "";
    }
    for (int i = 0; i < maxna5; i++) {
        ka5[i] = "";
    }
    for (int i = 0; i < maxna4; i++) {
        ka4[i] = "";
    }
    for (int i = 0; i < maxna3; i++) {
        ka3[i] = "";
    }
    for (int i = 0; i < maxnap; i++) {
        kap[i] = "";
    }
    for (int i = 0; i < maxnaf; i++) {
        kaf[i] = "";
    }
    for (int i = 0; i < maxnsb; i++) {
        ksb[i] = "";
    }
    for (int i = 0; i < maxnu; i++) {
        ku[i] = "";
    }
    for (int i = 0; i < maxnopb; i++) {
        kopb[i] = "";
    }
    for (int i = 0; i < maxnopd; i++) {
        kopd[i] = "";
    }
    for (int i = 0; i < maxndi; i++) {
        kdi[i] = "";
    }
    for (int i = 0; i < maxnti; i++) {
        kti[i] = "";
    }
    for (int i = 0; i < maxnt; i++) {
        kt[i] = "";
    }
    for (int i = 0; i < maxnt5; i++) {
        kt5[i] = "";
    }
    for (int i = 0; i < maxnt4; i++) {
        kt4[i] = "";
    }
    for (int i = 0; i < maxnpt; i++) {
        kpt[i] = "";
    }
    for (int i = 0; i < maxnbt; i++) {
        kbt[i] = "";
    }
    for (int i = 0; i < maxnat; i++) {
        kat[i] = "";
    }
    for (int i = 0; i < maxntt; i++) {
        ktt[i] = "";
    }
    for (int i = 0; i < maxnd; i++) {
        kd[i] = "";
    }
    for (int i = 0; i < maxnd5; i++) {
        kd5[i] = "";
    }
    for (int i = 0; i < maxnd4; i++) {
        kd4[i] = "";
    }
    for (int i = 0; i < maxnd3; i++) {
        kd3[i] = "";
    }
    for (int i = 0; i < maxnmp; i++) {
        kmp[i] = "";
    }
    for (int i = 0; i < maxnpp; i++) {
        kppr[i] = "";
    }
    for (int i = 0; i < maxncfb; i++) {
        kcfb[i] = "";
    }
    for (int i = 0; i < maxncfa; i++) {
        kcfa[i] = "";
    }
    for (int i = 0; i < maxnpi; i++) {
        kpi[i] = "";
    }
    for (int i = 0; i < maxnpi5; i++) {
        kpi5[i] = "";
    }
    for (int i = 0; i < maxnpi4; i++) {
        kpi4[i] = "";
    }

    // perform dynamic allocation of some global arrays
    atmcls.resize(maxtyp, -1);
    atmnum.resize(maxtyp, 0);
    ligand.resize(maxtyp);
    weight.resize(maxtyp);
    symbol.resize(maxtyp);
    describe.resize(maxtyp);
    anan.resize(maxclass, std::vector<double>(3));
    rad.resize(maxtyp);
    eps.resize(maxtyp);
    rad4.resize(maxtyp);
    eps4.resize(maxtyp);
    reduct.resize(maxtyp);
    prsiz.resize(maxclass);
    prdmp.resize(maxclass);
    prele.resize(maxclass);
    dspsix.resize(maxclass);
    dspdmp.resize(maxclass);
    chg.resize(maxtyp);
    cpele.resize(maxclass);
    cpalp.resize(maxclass);
    polr.resize(maxtyp);
    athl.resize(maxtyp);
    dthl.resize(maxtyp);
    pgrp.resize(maxtyp, std::vector<int>(maxval));
    pepk.resize(maxclass);
    peppre.resize(maxclass);
    pepdmp.resize(maxclass);
    pepl.resize(maxclass);
    ctchg.resize(maxclass);
    ctdmp.resize(maxclass);
    pbr.resize(maxtyp);
    csr.resize(maxtyp);
    gkr.resize(maxtyp);
    electron.resize(maxclass);
    ionize.resize(maxclass);
    repulse.resize(maxclass);
    biotyp.resize(maxbio);

    // initialize values of force field model parameters
    forcefield = "";
    for (int i = 0; i < maxtyp; i++) {
        ligand[i] = 0;
        weight[i] = 0.;
        symbol[i] = "";
        describe[i] = "";
        rad[i] = 0.;
        eps[i] = 0.;
        rad4[i] = 0.;
        eps4[i] = 0.;
        reduct[i] = 0.;
        chg[i] = 0.;
        polr[i] = 0.;
        athl[i] = 0.;
        dthl[i] = 0.;
        for (int j = 0; j < maxval; j++) {
            pgrp[i][j] = 0;
        }
        pbr[i] = 0.;
        csr[i] = 0.;
        gkr[i] = 0.;
    }
    for (int i = 0; i < maxclass; i++) {
        for (int j = 0; j < 3; j++) {
            anan[i][j] = 0.;
        }
        prsiz[i] = 0.;
        prdmp[i] = 0.;
        prele[i] = 0.;
        dspsix[i] = 0.;
        dspdmp[i] = 0.;
        cpele[i] = 0.;
        cpalp[i] = 0.;
        pepk[i] = 0.;
        peppre[i] = 0.;
        pepdmp[i] = 0.;
        pepl[i] = false;
        ctchg[i] = 0.;
        ctdmp[i] = 0.;
        electron[i] = 0.;
        ionize[i] = 0.;
        repulse[i] = 0.;
    }
    for (int i = 0; i < maxbio; i++) {
        biotyp[i] = 0;
    }

    // set default control parameters for local geometry terms
    bndtyp = "HARMONIC";
    bndunit = 1.;
    cbnd = 0.;
    qbnd = 0.;
    angunit = 1. / std::pow(radian,2);
    cang = 0.;
    qang = 0.;
    pang = 0.;
    sang = 0.;
    stbnunit = 1. / radian;
    ureyunit = 1.;
    cury = 0.;
    qury = 0.;
    aaunit = 1. / std::pow(radian,2);
    opbtyp = "W-D-C";
    opbunit = 1. / std::pow(radian,2);
    copb = 0.;
    qopb = 0.;
    popb = 0.;
    sopb = 0.;
    opdunit = 1.;
    copd = 0.;
    qopd = 0.;
    popd = 0.;
    sopd = 0.;
    idihunit = 1. / std::pow(radian,2);
    itorunit = 1.;
    torsunit = 1.;
    ptorunit = 1.;
    storunit = 1.;
    atorunit = 1. / radian;
    ttorunit = 1.;

    // set default control parameters for van der Waals terms
    vdwindex = "CLASS";
    vdwtyp = "LENNARD-JONES";
    radrule = "ARITHMETIC";
    radtyp = "R-MIN";
    radsiz = "RADIUS";
    epsrule = "GEOMETRIC";
    gausstyp = "NONE";
    ngauss = 0;
    abuck = 0.;
    bbuck = 0.;
    cbuck = 0.;
    ghal = 0.12;
    dhal = 0.07;
    v2scale = 0.;
    v3scale = 0.;
    v4scale = 1.;
    v5scale = 1.;
    use_vcorr = false;

    // set default control parameters for repulsion terms
    r2scale = 0.;
    r3scale = 0.;
    r4scale = 1.;
    r5scale = 1.;

    // set default control parameters for dispersion terms
    dsp2scale = 0.;
    dsp3scale = 0.;
    dsp4scale = 1.;
    dsp5scale = 1.;
    use_dcorr = false;

    // set default control parameters for charge-charge terms
    electric = coulomb;
    dielec = 1.;
    ebuffer = 0.;
    c1scale = 0.;
    c2scale = 0.;
    c3scale = 0.;
    c4scale = 1.;
    c5scale = 1.;
    neutnbr = false;
    neutcut = false;
    use_exfld = false;
    for (int i = 0; i < 3; i++) {
        exfld[i] = 0.;
    }

    // set default control parameters for atomic multipole terms
    pentyp = "GORDON1";
    m2scale = 0.;
    m3scale = 0.;
    m4scale = 1.;
    m5scale = 1.;
    d1scale = 0.;
    d2scale = 1.;
    d3scale = 1.;
    d4scale = 1.;
    use_chgpen = false;

    // set default control parameters for polarization terms
    poltyp = "MUTUAL";
    scrtyp = "S2U";
    politer = 100;
    poleps = 0.000001;
    uaccel = 2.;
    p2scale = 0.;
    p3scale = 0.;
    p4scale = 1.;
    p5scale = 1.;
    p2iscale = 0.;
    p3iscale = 0.;
    p4iscale = 0.5;
    p5iscale = 1.;
    u1scale = 1.;
    u2scale = 1.;
    u3scale = 1.;
    u4scale = 1.;
    w2scale = 1.;
    w3scale = 1.;
    w4scale = 1.;
    w5scale = 1.;
    use_thole = true;
    use_tholed = false;
    use_pred = false;
    use_ielscf = false;
    dpequal = false;
    use_expol = false;

    // set default control parameters for charge transfer terms
    ctrntyp = "SEPARATE";

    // set default control parameters for implicit solvation
    solvtyp = "";
    borntyp = "";

    // set default control parameters for reaction field
    rfsize = 1000000.;
    rfbulkd = 80.;
    rfterms = 1;

    // initialize some Merck Molecular force field parameters
    initmmff();
}


/////////////////////////////////////////////////////////
//                                                     //
//  initmmff.cpp  --  initialize some MMFF parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// "initmmff" initializes some parameter values for the Merck
// Molecular force field


#include "merck.h"

void initmmff()
{
    // perform dynamic allocation of some global arrays
    mmff_ka.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ka1.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ka2.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ka3.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ka4.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ka5.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ka6.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ka7.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ka8.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ang0.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ang1.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ang2.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ang3.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ang4.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ang5.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ang6.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ang7.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    mmff_ang8.resize(101, std::vector<std::vector<double>>(100, std::vector<double>(101,1000.)));
    stbn_abc.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_abc1.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba1.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_abc2.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba2.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_abc3.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba3.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_abc4.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba4.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_abc5.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba5.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_abc6.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba6.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_abc7.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba7.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_abc8.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba8.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_abc9.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba9.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_abc10.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba10.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_abc11.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));
    stbn_cba11.resize(100, std::vector<std::vector<double>>(100, std::vector<double>(100,1000.)));

    // initialize values for MMFF atom class equivalencies
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 500; j++) {
            eqclass[i][j] = 1000;
        }
    }

    // initialize values for MMFF aromatic ring parameters
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < maxtyp; j++) {
            mmffarom[i][j] = 0;
            mmffaromc[i][j] = 0;
            mmffaroma[i][j] = 0;
        }
    }

    // initialize values for MMFF bond stretching parameters
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            mmff_kb[i][j] = 1000.;
            mmff_kb1[i][j] = 1000.;
            mmff_b0[i][j] = 1000.;
            mmff_b1[i][j] = 1000.;
        }
    }

    // // initialize values for MMFF angle bending parameters
    // for (int i = 0; i < 101; i++) {
    //     for (int j = 0; j < 100; j++) {
    //         for (int k = 0; k < 101; k++) {
    //             mmff_ka[i][j][k] = 1000.;
    //             mmff_ka1[i][j][k] = 1000.;
    //             mmff_ka2[i][j][k] = 1000.;
    //             mmff_ka3[i][j][k] = 1000.;
    //             mmff_ka4[i][j][k] = 1000.;
    //             mmff_ka5[i][j][k] = 1000.;
    //             mmff_ka6[i][j][k] = 1000.;
    //             mmff_ka7[i][j][k] = 1000.;
    //             mmff_ka8[i][j][k] = 1000.;
    //             mmff_ang0[i][j][k] = 1000.;
    //             mmff_ang1[i][j][k] = 1000.;
    //             mmff_ang2[i][j][k] = 1000.;
    //             mmff_ang3[i][j][k] = 1000.;
    //             mmff_ang4[i][j][k] = 1000.;
    //             mmff_ang5[i][j][k] = 1000.;
    //             mmff_ang6[i][j][k] = 1000.;
    //             mmff_ang7[i][j][k] = 1000.;
    //             mmff_ang8[i][j][k] = 1000.;
    //         }
    //     }
    // }

    // // initialize values for MMFF stretch-bend parameters
    // for (int i = 0; i < 100; i++) {
    //     for (int j = 0; j < 100; j++) {
    //         for (int k = 0; k < 100; k++) {
    //             stbn_abc[i][j][k] = 1000.;
    //             stbn_cba[i][j][k] = 1000.;
    //             stbn_abc1[i][j][k] = 1000.;
    //             stbn_cba1[i][j][k] = 1000.;
    //             stbn_abc2[i][j][k] = 1000.;
    //             stbn_cba2[i][j][k] = 1000.;
    //             stbn_abc3[i][j][k] = 1000.;
    //             stbn_cba3[i][j][k] = 1000.;
    //             stbn_abc4[i][j][k] = 1000.;
    //             stbn_cba4[i][j][k] = 1000.;
    //             stbn_abc5[i][j][k] = 1000.;
    //             stbn_cba5[i][j][k] = 1000.;
    //             stbn_abc6[i][j][k] = 1000.;
    //             stbn_cba6[i][j][k] = 1000.;
    //             stbn_abc7[i][j][k] = 1000.;
    //             stbn_cba7[i][j][k] = 1000.;
    //             stbn_abc8[i][j][k] = 1000.;
    //             stbn_cba8[i][j][k] = 1000.;
    //             stbn_abc9[i][j][k] = 1000.;
    //             stbn_cba9[i][j][k] = 1000.;
    //             stbn_abc10[i][j][k] = 1000.;
    //             stbn_cba10[i][j][k] = 1000.;
    //             stbn_abc11[i][j][k] = 1000.;
    //             stbn_cba11[i][j][k] = 1000.;
    //         }
    //     }
    // }

    // initialize values for MMFF torsional parameters
    for (int i = 0; i < maxnt; i++) {
        kt[i] = "";
        kt_1[i] = "";
        kt_2[i] = "";
        t1[i][0] = 1000.;
        t1[i][1] = 1000.;
        t2[i][0] = 1000.;
        t2[i][1] = 1000.;
        t3[i][0] = 1000.;
        t3[i][1] = 1000.;
        t1_1[i][0] = 1000.;
        t1_1[i][1] = 1000.;
        t2_1[i][0] = 1000.;
        t2_1[i][1] = 1000.;
        t3_1[i][0] = 1000.;
        t3_1[i][1] = 1000.;
        t1_2[i][0] = 1000.;
        t1_2[i][1] = 1000.;
        t2_2[i][0] = 1000.;
        t2_2[i][1] = 1000.;
        t3_2[i][0] = 1000.;
        t3_2[i][1] = 1000.;
    }

    // initialize values for MMFF bond charge increment parameters
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            bci[i][j] = 1000.;
            bci_1[i][j] = 1000.;
        }
    }
}
