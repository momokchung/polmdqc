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
#include <string>

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
    atmcls.resize(maxtyp);
    atmnum.resize(maxtyp);
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
    pgrp.resize(maxval, std::vector<int>(maxtyp));
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
        atmcls[i] = 0;
        atmnum[i] = 0;
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
}

// c

// c

// c
// c     set default control parameters for local geometry terms
// c
//       bndtyp = 'HARMONIC'
//       bndunit = 1.0d0
//       cbnd = 0.0d0
//       qbnd = 0.0d0
//       angunit = 1.0d0 / radian**2
//       cang = 0.0d0
//       qang = 0.0d0
//       pang = 0.0d0
//       sang = 0.0d0
//       stbnunit = 1.0d0 / radian
//       ureyunit = 1.0d0
//       cury = 0.0d0
//       qury = 0.0d0
//       aaunit = 1.0d0 / radian**2
//       opbtyp = 'W-D-C'
//       opbunit = 1.0d0 / radian**2
//       copb = 0.0d0
//       qopb = 0.0d0
//       popb = 0.0d0
//       sopb = 0.0d0
//       opdunit = 1.0d0
//       copd = 0.0d0
//       qopd = 0.0d0
//       popd = 0.0d0
//       sopd = 0.0d0
//       idihunit = 1.0d0 / radian**2
//       itorunit = 1.0d0
//       torsunit = 1.0d0
//       ptorunit = 1.0d0
//       storunit = 1.0d0
//       atorunit = 1.0d0 / radian
//       ttorunit = 1.0d0
// c
// c     set default control parameters for van der Waals terms
// c
//       vdwindex = 'CLASS'
//       vdwtyp = 'LENNARD-JONES'
//       radrule = 'ARITHMETIC'
//       radtyp = 'R-MIN'
//       radsiz = 'RADIUS'
//       epsrule = 'GEOMETRIC'
//       gausstyp = 'NONE'
//       ngauss = 0
//       abuck = 0.0d0
//       bbuck = 0.0d0
//       cbuck = 0.0d0
//       ghal = 0.12d0
//       dhal = 0.07d0
//       v2scale = 0.0d0
//       v3scale = 0.0d0
//       v4scale = 1.0d0
//       v5scale = 1.0d0
//       use_vcorr = .false.
// c
// c     set default control parameters for repulsion terms
// c
//       r2scale = 0.0d0
//       r3scale = 0.0d0
//       r4scale = 1.0d0
//       r5scale = 1.0d0
// c
// c     set default control parameters for dispersion terms
// c
//       dsp2scale = 0.0d0
//       dsp3scale = 0.0d0
//       dsp4scale = 1.0d0
//       dsp5scale = 1.0d0
//       use_dcorr = .false.
// c
// c     set default control parameters for charge-charge terms
// c
//       electric = coulomb
//       dielec = 1.0d0
//       ebuffer = 0.0d0
//       c1scale = 0.0d0
//       c2scale = 0.0d0
//       c3scale = 0.0d0
//       c4scale = 1.0d0
//       c5scale = 1.0d0
//       neutnbr = .false.
//       neutcut = .false.
//       use_exfld = .false.
//       do i = 1, 3
//          exfld(i) = 0.0d0
//       end do
// c
// c     set default control parameters for atomic multipole terms
// c
//       pentyp = 'GORDON1'
//       m2scale = 0.0d0
//       m3scale = 0.0d0
//       m4scale = 1.0d0
//       m5scale = 1.0d0
//       d1scale = 0.0d0
//       d2scale = 1.0d0
//       d3scale = 1.0d0
//       d4scale = 1.0d0
//       use_chgpen = .false.
// c
// c     set default control parameters for polarization terms
// c
//       poltyp = 'MUTUAL'
//       scrtyp = 'S2U'
//       politer = 100
//       poleps = 0.000001d0
//       uaccel = 2.0d0
//       p2scale = 0.0d0
//       p3scale = 0.0d0
//       p4scale = 1.0d0
//       p5scale = 1.0d0
//       p2iscale = 0.0d0
//       p3iscale = 0.0d0
//       p4iscale = 0.5d0
//       p5iscale = 1.0d0
//       u1scale = 1.0d0
//       u2scale = 1.0d0
//       u3scale = 1.0d0
//       u4scale = 1.0d0
//       w2scale = 1.0d0
//       w3scale = 1.0d0
//       w4scale = 1.0d0
//       w5scale = 1.0d0
//       use_thole = .true.
//       use_tholed = .false.
//       use_pred = .false.
//       use_ielscf = .false.
//       dpequal = .false.
//       use_expol = .false.
// c
// c     set default control parameters for charge transfer terms
// c
//       ctrntyp = 'SEPARATE'
// c
// c     set default control parameters for implicit solvation
// c
//       solvtyp = ""
//       borntyp = ""
// c
// c     set default control parameters for reaction field
// c
//       rfsize = 1000000.0d0
//       rfbulkd = 80.0d0
//       rfterms = 1
// c
// c     initialize some Merck Molecular force field parameters
// c
//       call initmmff
//       return
//       end
// c
// c
// c     ################################################################
// c     ##                                                            ##
// c     ##  subroutine initmmff  --  initialize some MMFF parameters  ##
// c     ##                                                            ##
// c     ################################################################
// c
// c
// c     "initmmff" initializes some parameter values for the Merck
// c     Molecular force field
// c
// c
//       subroutine initmmff
//       use ktorsn
//       use merck
//       implicit none
//       integer i,j,k
//       character*16 blank16
// c
// c
// c     define blank character strings of various lengths
// c
//       blank16 = '                '
// c
// c     perform dynamic allocation of some global arrays
// c
//       if (.not. allocated(mmff_ka))
//      &   allocate (mmff_ka(0:100,100,0:100))
//       if (.not. allocated(mmff_ka1))
//      &   allocate (mmff_ka1(0:100,100,0:100))
//       if (.not. allocated(mmff_ka2))
//      &   allocate (mmff_ka2(0:100,100,0:100))
//       if (.not. allocated(mmff_ka3))
//      &   allocate (mmff_ka3(0:100,100,0:100))
//       if (.not. allocated(mmff_ka4))
//      &   allocate (mmff_ka4(0:100,100,0:100))
//       if (.not. allocated(mmff_ka5))
//      &   allocate (mmff_ka5(0:100,100,0:100))
//       if (.not. allocated(mmff_ka6))
//      &   allocate (mmff_ka6(0:100,100,0:100))
//       if (.not. allocated(mmff_ka7))
//      &   allocate (mmff_ka7(0:100,100,0:100))
//       if (.not. allocated(mmff_ka8))
//      &   allocate (mmff_ka8(0:100,100,0:100))
//       if (.not. allocated(mmff_ang0))
//      &   allocate (mmff_ang0(0:100,100,0:100))
//       if (.not. allocated(mmff_ang1))
//      &   allocate (mmff_ang1(0:100,100,0:100))
//       if (.not. allocated(mmff_ang2))
//      &   allocate (mmff_ang2(0:100,100,0:100))
//       if (.not. allocated(mmff_ang3))
//      &   allocate (mmff_ang3(0:100,100,0:100))
//       if (.not. allocated(mmff_ang4))
//      &   allocate (mmff_ang4(0:100,100,0:100))
//       if (.not. allocated(mmff_ang5))
//      &   allocate (mmff_ang5(0:100,100,0:100))
//       if (.not. allocated(mmff_ang6))
//      &   allocate (mmff_ang6(0:100,100,0:100))
//       if (.not. allocated(mmff_ang7))
//      &   allocate (mmff_ang7(0:100,100,0:100))
//       if (.not. allocated(mmff_ang8))
//      &   allocate (mmff_ang8(0:100,100,0:100))
//       if (.not. allocated(stbn_abc))
//      &   allocate (stbn_abc(100,100,100))
//       if (.not. allocated(stbn_cba))
//      &   allocate (stbn_cba(100,100,100))
//       if (.not. allocated(stbn_abc1))
//      &   allocate (stbn_abc1(100,100,100))
//       if (.not. allocated(stbn_cba1))
//      &   allocate (stbn_cba1(100,100,100))
//       if (.not. allocated(stbn_abc2))
//      &   allocate (stbn_abc2(100,100,100))
//       if (.not. allocated(stbn_cba2))
//      &   allocate (stbn_cba2(100,100,100))
//       if (.not. allocated(stbn_abc3))
//      &   allocate (stbn_abc3(100,100,100))
//       if (.not. allocated(stbn_cba3))
//      &   allocate (stbn_cba3(100,100,100))
//       if (.not. allocated(stbn_abc4))
//      &   allocate (stbn_abc4(100,100,100))
//       if (.not. allocated(stbn_cba4))
//      &   allocate (stbn_cba4(100,100,100))
//       if (.not. allocated(stbn_abc5))
//      &   allocate (stbn_abc5(100,100,100))
//       if (.not. allocated(stbn_cba5))
//      &   allocate (stbn_cba5(100,100,100))
//       if (.not. allocated(stbn_abc6))
//      &   allocate (stbn_abc6(100,100,100))
//       if (.not. allocated(stbn_cba6))
//      &   allocate (stbn_cba6(100,100,100))
//       if (.not. allocated(stbn_abc7))
//      &   allocate (stbn_abc7(100,100,100))
//       if (.not. allocated(stbn_cba7))
//      &   allocate (stbn_cba7(100,100,100))
//       if (.not. allocated(stbn_abc8))
//      &   allocate (stbn_abc8(100,100,100))
//       if (.not. allocated(stbn_cba8))
//      &   allocate (stbn_cba8(100,100,100))
//       if (.not. allocated(stbn_abc9))
//      &   allocate (stbn_abc9(100,100,100))
//       if (.not. allocated(stbn_cba9))
//      &   allocate (stbn_cba9(100,100,100))
//       if (.not. allocated(stbn_abc10))
//      &   allocate (stbn_abc10(100,100,100))
//       if (.not. allocated(stbn_cba10))
//      &   allocate (stbn_cba10(100,100,100))
//       if (.not. allocated(stbn_abc11))
//      &   allocate (stbn_abc11(100,100,100))
//       if (.not. allocated(stbn_cba11))
//      &   allocate (stbn_cba11(100,100,100))
// c
// c     initialize values for MMFF atom class equivalencies
// c
//       do i = 1, 5
//          do j = 1, 500
//             eqclass(j,i) = 1000
//          end do
//       end do
// c
// c     initialize values for MMFF aromatic ring parameters
// c
//       do i = 1, 6
//          do j = 1, maxtyp
//             mmffarom(j,i) = 0
//             mmffaromc(j,i) = 0
//             mmffaroma(j,i) = 0
//          end do
//       end do
// c
// c     initialize values for MMFF bond stretching parameters
// c
//       do i = 1, 100
//          do j = 1, 100
//             mmff_kb(j,i) = 1000.0d0
//             mmff_kb1(j,i) = 1000.0d0
//             mmff_b0(j,i) = 1000.0d0
//             mmff_b1(j,i) = 1000.0d0
//          end do
//       end do
// c
// c     initialize values for MMFF angle bending parameters
// c
//       do i = 0, 100
//          do j = 1, 100
//             do k = 0, 100
//                mmff_ka(k,j,i) = 1000.0d0
//                mmff_ka1(k,j,i) = 1000.0d0
//                mmff_ka2(k,j,i) = 1000.0d0
//                mmff_ka3(k,j,i) = 1000.0d0
//                mmff_ka4(k,j,i) = 1000.0d0
//                mmff_ka5(k,j,i) = 1000.0d0
//                mmff_ka6(k,j,i) = 1000.0d0
//                mmff_ka7(k,j,i) = 1000.0d0
//                mmff_ka8(k,j,i) = 1000.0d0
//                mmff_ang0(k,j,i) = 1000.0d0
//                mmff_ang1(k,j,i) = 1000.0d0
//                mmff_ang2(k,j,i) = 1000.0d0
//                mmff_ang3(k,j,i) = 1000.0d0
//                mmff_ang4(k,j,i) = 1000.0d0
//                mmff_ang5(k,j,i) = 1000.0d0
//                mmff_ang6(k,j,i) = 1000.0d0
//                mmff_ang7(k,j,i) = 1000.0d0
//                mmff_ang8(k,j,i) = 1000.0d0
//             end do
//          end do
//       end do
// c
// c     initialize values for MMFF stretch-bend parameters
// c
//       do i = 1, 100
//          do j = 1, 100
//             do k = 1, 100
//                stbn_abc(k,j,i) = 1000.0d0
//                stbn_cba(k,j,i) = 1000.0d0
//                stbn_abc1(k,j,i) = 1000.0d0
//                stbn_cba1(k,j,i) = 1000.0d0
//                stbn_abc2(k,j,i) = 1000.0d0
//                stbn_cba2(k,j,i) = 1000.0d0
//                stbn_abc3(k,j,i) = 1000.0d0
//                stbn_cba3(k,j,i) = 1000.0d0
//                stbn_abc4(k,j,i) = 1000.0d0
//                stbn_cba4(k,j,i) = 1000.0d0
//                stbn_abc5(k,j,i) = 1000.0d0
//                stbn_cba5(k,j,i) = 1000.0d0
//                stbn_abc6(k,j,i) = 1000.0d0
//                stbn_cba6(k,j,i) = 1000.0d0
//                stbn_abc7(k,j,i) = 1000.0d0
//                stbn_cba7(k,j,i) = 1000.0d0
//                stbn_abc8(k,j,i) = 1000.0d0
//                stbn_cba8(k,j,i) = 1000.0d0
//                stbn_abc9(k,j,i) = 1000.0d0
//                stbn_cba9(k,j,i) = 1000.0d0
//                stbn_abc10(k,j,i) = 1000.0d0
//                stbn_cba10(k,j,i) = 1000.0d0
//                stbn_abc11(k,j,i) = 1000.0d0
//                stbn_cba11(k,j,i) = 1000.0d0
//             end do
//          end do
//       end do
// c
// c     initialize values for MMFF torsional parameters
// c
//       do i = 1, maxnt
//          kt(i) = blank16
//          kt_1(i) = blank16
//          kt_2(i) = blank16
//          t1(1,i) = 1000.0d0
//          t1(2,i) = 1000.0d0
//          t2(1,i) = 1000.0d0
//          t2(2,i) = 1000.0d0
//          t3(1,i) = 1000.0d0
//          t3(2,i) = 1000.0d0
//          t1_1(1,i) = 1000.0d0
//          t1_1(2,i) = 1000.0d0
//          t2_1(1,i) = 1000.0d0
//          t2_1(2,i) = 1000.0d0
//          t3_1(1,i) = 1000.0d0
//          t3_1(2,i) = 1000.0d0
//          t1_2(1,i) = 1000.0d0
//          t1_2(2,i) = 1000.0d0
//          t2_2(1,i) = 1000.0d0
//          t2_2(2,i) = 1000.0d0
//          t3_2(1,i) = 1000.0d0
//          t3_2(2,i) = 1000.0d0
//       end do
// c
// c     initialize values for MMFF bond charge increment parameters
// c
//       do i = 1, 100
//          do j = 1, 100
//             bci(j,i) = 1000.0d0
//             bci_1(j,i) = 1000.0d0
//          end do
//       end do
//       return
//       end
