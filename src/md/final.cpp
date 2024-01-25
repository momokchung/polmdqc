// Author: Moses KJ Chung
// Year:   2023

#include "darray.h"
#include "final.h"
#include "mod.h"

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  final  --  final actions before program exit  //
//                                                //
////////////////////////////////////////////////////

// "final" performs any final program actions such as deallocation
// of global memory, prints a status message, and then pauses if
// necessary to avoid closing the execution window

void final()
{
    // print a final status message before exiting PolMDQC
    if (debug) {
        printf("\n PolMDQC is Exiting following Normal Termination\n");
    }

    // deallocation of global arrays from module align
    ifit.deallocate();
    wfit.deallocate();

    // deallocation of global arrays from module alphmol
    surf.deallocate();
    vol.deallocate();
    mean.deallocate();
    gauss.deallocate();
    dsurf.deallocate();
    dvol.deallocate();
    dmean.deallocate();
    dgauss.deallocate();
    ndsurf.deallocate();
    ndvol.deallocate();
    ndmean.deallocate();
    ndgauss.deallocate();

    // deallocation of global arrays from module analyz
    aesum.deallocate();
    aeb.deallocate();
    aea.deallocate();
    aeba.deallocate();
    aeub.deallocate();
    aeaa.deallocate();
    aeopb.deallocate();
    aeopd.deallocate();
    aeid.deallocate();
    aeit.deallocate();
    aet.deallocate();
    aept.deallocate();
    aebt.deallocate();
    aeat.deallocate();
    aett.deallocate();
    aev.deallocate();
    aer.deallocate();
    aedsp.deallocate();
    aec.deallocate();
    aecd.deallocate();
    aed.deallocate();
    aem.deallocate();
    aep.deallocate();
    aect.deallocate();
    aerxf.deallocate();
    aes.deallocate();
    aelf.deallocate();
    aeg.deallocate();
    aex.deallocate();

    // deallocation of global arrays from module angbnd
    iang.deallocate();
    ak.deallocate();
    anat.deallocate();
    afld.deallocate();

    // deallocation of global arrays from module angpot
    angtyp.deallocate();

    // deallocation of global arrays from module atmlst
    bndlist.deallocate();
    anglist.deallocate();
    balist.deallocate();

    // deallocation of global arrays from module atomid
    tag.deallocate();
    atomClass.deallocate();
    atomic.deallocate();
    valence.deallocate();
    mass.deallocate();
    name.deallocate();
    story.deallocate();

    // deallocation of global arrays from module atoms
    type.deallocate();
    x.deallocate();
    y.deallocate();
    z.deallocate();

    // deallocation of global arrays from module bitor
    ibitor.deallocate();

    // deallocation of global arrays from module bndstr
    ibnd.deallocate();
    bk.deallocate();
    bl.deallocate();

    // deallocation of global arrays from module cell
    icell.deallocate();

    // deallocation of global arrays from module chgpen
    pcore.deallocate();
    pval.deallocate();
    pval0.deallocate();
    palpha.deallocate();

    // deallocation of global arrays from module chgtrn
    chgct.deallocate();
    dmpct.deallocate();

    // deallocation of global arrays from module chunks
    pmetable.deallocate();

    // deallocation of global arrays from module couple
    n12.deallocate();
    n13.deallocate();
    n14.deallocate();
    n15.deallocate();
    i12.deallocate();
    i13.deallocate();
    i14.deallocate();
    i15.deallocate();

    // deallocation of global arrays from module deriv
    desum.deallocate();
    deb.deallocate();
    dea.deallocate();
    deba.deallocate();
    deub.deallocate();
    deaa.deallocate();
    deopb.deallocate();
    deopd.deallocate();
    deid.deallocate();
    deit.deallocate();
    det.deallocate();
    dept.deallocate();
    debt.deallocate();
    deat.deallocate();
    dett.deallocate();
    dev.deallocate();
    der.deallocate();
    dedsp.deallocate();
    dec.deallocate();
    decd.deallocate();
    ded.deallocate();
    dem.deallocate();
    dep.deallocate();
    dect.deallocate();
    derxf.deallocate();
    des.deallocate();
    delf.deallocate();
    deg.deallocate();
    dex.deallocate();
    te.deallocate();
    ndesum.deallocate();
    ndeb.deallocate();
    ndea.deallocate();
    ndeba.deallocate();
    ndeub.deallocate();
    ndeaa.deallocate();
    ndeopb.deallocate();
    ndeopd.deallocate();
    ndeid.deallocate();
    ndeit.deallocate();
    ndet.deallocate();
    ndept.deallocate();
    ndebt.deallocate();
    ndeat.deallocate();
    ndett.deallocate();
    ndev.deallocate();
    nder.deallocate();
    ndedsp.deallocate();
    ndec.deallocate();
    ndecd.deallocate();
    nded.deallocate();
    ndem.deallocate();
    ndep.deallocate();
    ndect.deallocate();
    nderxf.deallocate();
    ndes.deallocate();
    ndelf.deallocate();
    ndeg.deallocate();
    ndex.deallocate();

    // deallocation of global arrays from module energi
    escale.deallocate();

    // deallocation of global arrays from module expol
    kpep.deallocate();
    prepep.deallocate();
    dmppep.deallocate();
    polscale.deallocate();
    polinv.deallocate();
    lpep.deallocate();

    // deallocation of global arrays from module fft
    ffttable.deallocate();

    // deallocation of global arrays from module fields
    biotyp.deallocate();

    // deallocation of global arrays from module group
    kgrp.deallocate();
    grplist.deallocate();
    igrp.deallocate();
    grpmass.deallocate();
    wgrp.deallocate();

    // deallocation of global arrays from module ielscf
    uaux.deallocate();
    upaux.deallocate();
    vaux.deallocate();
    vpaux.deallocate();
    aaux.deallocate();
    apaux.deallocate();

    // deallocation of global arrays from module kanang
    anan.deallocate();

    // deallocation of global arrays from module kangs
    acon.deallocate();
    acon5.deallocate();
    acon4.deallocate();
    acon3.deallocate();
    aconp.deallocate();
    aconf.deallocate();
    ang.deallocate();
    ang5.deallocate();
    ang4.deallocate();
    ang3.deallocate();
    angp.deallocate();
    angf.deallocate();
    ka.deallocate();
    ka5.deallocate();
    ka4.deallocate();
    ka3.deallocate();
    kap.deallocate();
    kaf.deallocate();

    // deallocation of global arrays from module kantor
    atcon.deallocate();
    kat.deallocate();

    // deallocation of global arrays from module katoms
    atmcls.deallocate();
    atmnum.deallocate();
    ligand.deallocate();
    weight.deallocate();
    symbol.deallocate();
    describe.deallocate();

    // deallocation of global arrays from module kbonds
    bcon.deallocate();
    bcon5.deallocate();
    bcon4.deallocate();
    bcon3.deallocate();
    blen.deallocate();
    blen5.deallocate();
    blen4.deallocate();
    blen3.deallocate();
    dlen.deallocate();
    kb.deallocate();
    kb5.deallocate();
    kb4.deallocate();
    kb3.deallocate();
    kel.deallocate();

    // deallocation of global arrays from module kcflux
    cflb.deallocate();
    cfla.deallocate();
    cflab.deallocate();
    kcfb.deallocate();
    kcfa.deallocate();

    // deallocation of global arrays from module kchrge
    chg.deallocate();

    // deallocation of global arrays from module kcpen
    cpele.deallocate();
    cpalp.deallocate();

    // deallocation of global arrays from module kctrn
    ctchg.deallocate();
    ctdmp.deallocate();

    // deallocation of global arrays from module kdipol
    dpl.deallocate();
    dpl5.deallocate();
    dpl4.deallocate();
    dpl3.deallocate();
    pos.deallocate();
    pos5.deallocate();
    pos4.deallocate();
    pos3.deallocate();
    kd.deallocate();
    kd5.deallocate();
    kd4.deallocate();
    kd3.deallocate();

    // deallocation of global arrays from module kdsp
    dspsix.deallocate();
    dspdmp.deallocate();

    // deallocation of global arrays from module kexpl
    pepk.deallocate();
    peppre.deallocate();
    pepdmp.deallocate();
    pepl.deallocate();

    // deallocation of global arrays from module khbond
    radhb.deallocate();
    epshb.deallocate();
    khb.deallocate();

    // deallocation of global arrays from module kiprop
    dcon.deallocate();
    tdi.deallocate();
    kdi.deallocate();

    // deallocation of global arrays from module kitors
    ti1.deallocate();
    ti2.deallocate();
    ti3.deallocate();
    kti.deallocate();

    // deallocation of global arrays from module kmulti
    multip.deallocate();
    mpaxis.deallocate();
    kmp.deallocate();

    // deallocation of global arrays from module kopbnd
    opbn.deallocate();
    kopb.deallocate();

    // deallocation of global arrays from module kopdst
    opds.deallocate();
    kopd.deallocate();

    // deallocation of global arrays from module korbs
    electron.deallocate();
    ionize.deallocate();
    repulse.deallocate();
    sslope.deallocate();
    sslope5.deallocate();
    sslope4.deallocate();
    tslope.deallocate();
    tslope5.deallocate();
    tslope4.deallocate();
    kpi.deallocate();
    kpi5.deallocate();
    kpi4.deallocate();

    // deallocation of global arrays from module kpitor
    ptcon.deallocate();
    kpt.deallocate();

    // deallocation of global arrays from module kpolpr
    thlpr.deallocate();
    thdpr.deallocate();
    kppr.deallocate();

    // deallocation of global arrays from module kpolr
    polr.deallocate();
    athl.deallocate();
    dthl.deallocate();
    pgrp.deallocate();

    // deallocation of global arrays from module krepl
    prsiz.deallocate();
    prdmp.deallocate();
    prele.deallocate();

    // deallocation of global arrays from module ksolut
    pbr.deallocate();
    csr.deallocate();
    gkr.deallocate();

    // deallocation of global arrays from module kstbnd
    stbn.deallocate();
    ksb.deallocate();

    // deallocation of global arrays from module ksttor
    btcon.deallocate();
    kbt.deallocate();

    // deallocation of global arrays from module ktorsn
    t1.deallocate();
    t2.deallocate();
    t3.deallocate();
    t4.deallocate();
    t5.deallocate();
    t6.deallocate();
    t15.deallocate();
    t25.deallocate();
    t35.deallocate();
    t45.deallocate();
    t55.deallocate();
    t65.deallocate();
    t14.deallocate();
    t24.deallocate();
    t34.deallocate();
    t44.deallocate();
    t54.deallocate();
    t64.deallocate();
    kt.deallocate();
    kt5.deallocate();
    kt4.deallocate();

    // deallocation of global arrays from module ktrtor
    tnx.deallocate();
    tny.deallocate();
    ttx.deallocate();
    tty.deallocate();
    tbf.deallocate();
    tbx.deallocate();
    tby.deallocate();
    tbxy.deallocate();
    ktt.deallocate();

    // deallocation of global arrays from module kurybr
    ucon.deallocate();
    dst13.deallocate();
    ku.deallocate();

    // deallocation of global arrays from module kvdwpr
    radpr.deallocate();
    epspr.deallocate();
    kvpr.deallocate();

    // deallocation of global arrays from module kvdws
    rad.deallocate();
    eps.deallocate();
    rad4.deallocate();
    eps4.deallocate();
    reduct.deallocate();

    // deallocation of global arrays from module molcul
    imol.deallocate();
    kmol.deallocate();
    molcule.deallocate();
    molmass.deallocate();

    // deallocation of global arrays from module mpole
    ipole.deallocate();
    polsiz.deallocate();
    pollist.deallocate();
    zaxis.deallocate();
    xaxis.deallocate();
    yaxis.deallocate();
    mono0.deallocate();
    pole.deallocate();
    rpole.deallocate();
    polaxe.deallocate();

    // deallocation of global arrays from module mutant
    imut.deallocate();
    type0.deallocate();
    class0.deallocate();
    type1.deallocate();
    class1.deallocate();
    mut.deallocate();

    // deallocation of global arrays from module neigh
    nvlst.deallocate();
    vlst.deallocate();
    nelst.deallocate();
    elst.deallocate();
    nulst.deallocate();
    ulst.deallocate();
    xvold.deallocate();
    yvold.deallocate();
    zvold.deallocate();
    xeold.deallocate();
    yeold.deallocate();
    zeold.deallocate();
    xuold.deallocate();
    yuold.deallocate();
    zuold.deallocate();

    // deallocation of global arrays from module nonpol
    radcav.deallocate();
    raddsp.deallocate();
    epsdsp.deallocate();
    cdsp.deallocate();

    // deallocation of global arrays from module pdb
    resnum.deallocate();
    resatm.deallocate();
    npdb12.deallocate();
    ipdb12.deallocate();
    pdblist.deallocate();
    xpdb.deallocate();
    ypdb.deallocate();
    zpdb.deallocate();
    pdbres.deallocate();
    pdbsym.deallocate();
    pdbatm.deallocate();
    pdbtyp.deallocate();

    // deallocation of global arrays from module pme
    igrid.deallocate();
    bsmod1.deallocate();
    bsmod2.deallocate();
    bsmod3.deallocate();
    bsbuild.deallocate();
    thetai1.deallocate();
    thetai2.deallocate();
    thetai3.deallocate();
    qgrid.deallocate();
    qfac.deallocate();

    // deallocation of global arrays from module polar
    ipolar.deallocate();
    jpolar.deallocate();
    polarity.deallocate();
    thole.deallocate();
    tholed.deallocate();
    pdamp.deallocate();
    thlval.deallocate();
    thdval.deallocate();
    udir.deallocate();
    udirp.deallocate();
    udirs.deallocate();
    udirps.deallocate();
    uind.deallocate();
    uinp.deallocate();
    uinds.deallocate();
    uinps.deallocate();
    uexact.deallocate();
    douind.deallocate();

    // deallocation of global arrays from module polgrp
    np11.deallocate();
    np12.deallocate();
    np13.deallocate();
    np14.deallocate();
    ip11.deallocate();
    ip12.deallocate();
    ip13.deallocate();
    ip14.deallocate();

    // deallocation of global arrays from module polopt
    copt.deallocate();
    copm.deallocate();
    uopt.deallocate();
    uoptp.deallocate();
    uopts.deallocate();
    uoptps.deallocate();
    fopt.deallocate();
    foptp.deallocate();

    // deallocation of global arrays from module polpcg
    mindex.deallocate();
    minv.deallocate();

    // deallocation of global arrays from module poltcg
    uad.deallocate();
    uap.deallocate();
    ubd.deallocate();
    ubp.deallocate();

    // deallocation of global arrays from module repel
    irep.deallocate();
    replist.deallocate();
    sizpr.deallocate();
    dmppr.deallocate();
    elepr.deallocate();
    repole.deallocate();
    rrepole.deallocate();

    // deallocation of global arrays from module restrn
    ipfix.deallocate();
    kpfix.deallocate();
    idfix.deallocate();
    iafix.deallocate();
    itfix.deallocate();
    igfix.deallocate();
    ichir.deallocate();
    xpfix.deallocate();
    ypfix.deallocate();
    zpfix.deallocate();
    pfix.deallocate();
    dfix.deallocate();
    afix.deallocate();
    tfix.deallocate();
    gfix.deallocate();
    chir.deallocate();

    // deallocation of global arrays from module rigid
    xrb.deallocate();
    yrb.deallocate();
    zrb.deallocate();
    rbc.deallocate();

    // deallocation of global arrays from module scales
    scale.deallocate();

    // deallocation of global arrays from module solute
    rsolv.deallocate();
    rdescr.deallocate();
    asolv.deallocate();
    rborn.deallocate();
    drb.deallocate();
    drbp.deallocate();
    drobc.deallocate();
    gpol.deallocate();
    shct.deallocate();
    aobc.deallocate();
    bobc.deallocate();
    gobc.deallocate();
    vsolv.deallocate();
    wace.deallocate();
    s2ace.deallocate();
    uace.deallocate();

    // deallocation of global arrays from module tarray
    tindex.deallocate();
    tdipdip.deallocate();

    // deallocation of global arrays from module tors
    itors.deallocate();
    tors1.deallocate();
    tors2.deallocate();
    tors3.deallocate();
    tors4.deallocate();
    tors5.deallocate();
    tors6.deallocate();

    // deallocation of global arrays from module uprior
    udalt.deallocate();
    upalt.deallocate();
    usalt.deallocate();
    upsalt.deallocate();

    // deallocation of global arrays from module usage
    iuse.deallocate();
    use.deallocate();

    // deallocation of global arrays from module vdw
    ivdw.deallocate();
    jvdw.deallocate();
    mvdw.deallocate();
    ired.deallocate();
    kred.deallocate();
    xred.deallocate();
    yred.deallocate();
    zred.deallocate();
    radmin.deallocate();
    epsilon.deallocate();
    radmin4.deallocate();
    epsilon4.deallocate();
    radhbnd.deallocate();
    epshbnd.deallocate();

    // deallocation of global arrays from module warp
    m2.deallocate();
}
}
