// Author: Moses KJ Chung
// Year:   2024

#include "alterchg.h"
#include "analyz.h"
#include "atoms.h"
#include "calcMode.h"
#include "deriv.h"
#include "empole.h"
#include "energi.h"
#include "energy.h"
#include "fatal.h"
#include "inter.h"
#include "openmp.h"
#include "potent.h"
#include "virial.h"
#include <cmath>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  energy  --  evaluates energy, gradient, and virial  //
//                                                      //
//////////////////////////////////////////////////////////

// "energy" calls the subroutines to calculate the potential
// energy, gradient and virial terms and sums their totals

template <CalcMode CalculationMode>
void energy()
{
    // choose calculation mode
    constexpr CalcFlag flags = getCalculationFlags<CalculationMode>();
    constexpr bool do_e = flags.do_energy;
    constexpr bool do_a = flags.do_analysis;
    constexpr bool do_g = flags.do_gradient;
    constexpr bool do_v = flags.do_virial;

    // zero out each of the potential energy components
    if constexpr (do_e) {
        esum = 0.;
        eb = 0.;
        ea = 0.;
        eba = 0.;
        eub = 0.;
        eaa = 0.;
        eopb = 0.;
        eopd = 0.;
        eid = 0.;
        eit = 0.;
        et = 0.;
        ept = 0.;
        ebt = 0.;
        eat = 0.;
        ett = 0.;
        ev = 0.;
        er = 0.;
        edsp = 0.;
        ec = 0.;
        ecd = 0.;
        ed = 0.;
        em = 0.;
        ep = 0.;
        ect = 0.;
        erxf = 0.;
        es = 0.;
        elf = 0.;
        eg = 0.;
        ex = 0.;
    }

    // zero out the total intermolecular energy
    if constexpr (do_a) einter = 0.;

    // perform dynamic allocation of escale
    escale.allocate(nthread,n);

    if constexpr (do_a) {
        // perform dynamic allocation of some global arrays
        aesum.allocate(n);
        aeb.allocate(n);
        aea.allocate(n);
        aeba.allocate(n);
        aeub.allocate(n);
        aeaa.allocate(n);
        aeopb.allocate(n);
        aeopd.allocate(n);
        aeid.allocate(n);
        aeit.allocate(n);
        aet.allocate(n);
        aept.allocate(n);
        aebt.allocate(n);
        aeat.allocate(n);
        aett.allocate(n);
        aev.allocate(n);
        aer.allocate(n);
        aedsp.allocate(n);
        aec.allocate(n);
        aecd.allocate(n);
        aed.allocate(n);
        aem.allocate(n);
        aep.allocate(n);
        aect.allocate(n);
        aerxf.allocate(n);
        aes.allocate(n);
        aelf.allocate(n);
        aeg.allocate(n);
        aex.allocate(n);

        // zero out energy partitioning components for each atom
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            aesum[i] = 0.;
            aeb[i] = 0.;
            aea[i] = 0.;
            aeba[i] = 0.;
            aeub[i] = 0.;
            aeaa[i] = 0.;
            aeopb[i] = 0.;
            aeopd[i] = 0.;
            aeid[i] = 0.;
            aeit[i] = 0.;
            aet[i] = 0.;
            aept[i] = 0.;
            aebt[i] = 0.;
            aeat[i] = 0.;
            aett[i] = 0.;
            aev[i] = 0.;
            aer[i] = 0.;
            aedsp[i] = 0.;
            aec[i] = 0.;
            aecd[i] = 0.;
            aed[i] = 0.;
            aem[i] = 0.;
            aep[i] = 0.;
            aect[i] = 0.;
            aerxf[i] = 0.;
            aes[i] = 0.;
            aelf[i] = 0.;
            aeg[i] = 0.;
            aex[i] = 0.;
        }
    }

    if constexpr (do_g) {
        // perform dynamic allocation of some global arrays
        desum.allocate(3*n);
        deb.allocate(3*n);
        dea.allocate(3*n);
        deba.allocate(3*n);
        deub.allocate(3*n);
        deaa.allocate(3*n);
        deopb.allocate(3*n);
        deopd.allocate(3*n);
        deid.allocate(3*n);
        deit.allocate(3*n);
        det.allocate(3*n);
        dept.allocate(3*n);
        debt.allocate(3*n);
        deat.allocate(3*n);
        dett.allocate(3*n);
        dev.allocate(3*n);
        der.allocate(3*n);
        dedsp.allocate(3*n);
        dec.allocate(3*n);
        decd.allocate(3*n);
        ded.allocate(3*n);
        dem.allocate(3*n);
        dep.allocate(3*n);
        dect.allocate(3*n);
        derxf.allocate(3*n);
        des.allocate(3*n);
        delf.allocate(3*n);
        deg.allocate(3*n);
        dex.allocate(3*n);
        te.allocate(3*n);

        // zero out each of the first derivative components
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                desum[3*i+j] = 0.;
                deb[3*i+j] = 0.;
                dea[3*i+j] = 0.;
                deba[3*i+j] = 0.;
                deub[3*i+j] = 0.;
                deaa[3*i+j] = 0.;
                deopb[3*i+j] = 0.;
                deopd[3*i+j] = 0.;
                deid[3*i+j] = 0.;
                deit[3*i+j] = 0.;
                det[3*i+j] = 0.;
                dept[3*i+j] = 0.;
                debt[3*i+j] = 0.;
                deat[3*i+j] = 0.;
                dett[3*i+j] = 0.;
                dev[3*i+j] = 0.;
                der[3*i+j] = 0.;
                dedsp[3*i+j] = 0.;
                dec[3*i+j] = 0.;
                decd[3*i+j] = 0.;
                ded[3*i+j] = 0.;
                dem[3*i+j] = 0.;
                dep[3*i+j] = 0.;
                dect[3*i+j] = 0.;
                derxf[3*i+j] = 0.;
                des[3*i+j] = 0.;
                delf[3*i+j] = 0.;
                deg[3*i+j] = 0.;
                dex[3*i+j] = 0.;
            }
        }
    }

    if constexpr (do_g and do_v) {
        // zero out the virial and the intermolecular energy
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                vir[i][j] = 0.;
            }
        }
    }

    // update the pairwise interaction neighbor lists
    // if (use_list) nblist();

    // remove any previous use of the replicates method
    // cutoff = 0.;
    // replica(cutoff);

    // remove any previous use of the replicates method
    // cutoff = 0.;
    // replica(cutoff);

    // many implicit solvation models require Born radii
    // if (use_born) born();

    // alter partial charges and multipoles for charge flux
    if (use_chgflx) alterchg();

    // modify bond and torsion constants for pisystem
    // if (use_orbit) picalc;

    // call the local geometry energy component routines
    // if (use_bond) ebond();
    // if (use_angle) eangle();
    // if (use_strbnd) estrbnd();
    // if (use_urey) eurey();
    // if (use_angang) eangang();
    // if (use_opbend) eopbend();
    // if (use_opdist) eopdist();
    // if (use_improp) eimprop();
    // if (use_imptor) eimptor();
    // if (use_tors) etors();
    // if (use_pitors) epitors();
    // if (use_strtor) estrtor();
    // if (use_angtor) eangtor();
    // if (use_tortor) etortor();

    // call the electrostatic energy component routines
    // if (use_charge) echarge();
    // if (use_chgdpl) echgdpl();
    // if (use_dipole) edipole();
    if (use_mpole) empole<CalculationMode>();
    // if (use_polar) epolar();
    // if (use_chgtrn) echgtrn();
    // if (use_rxnfld) erxnfld();

    // call the van der Waals energy component routines
    // if (use_vdw) {
    //     if (vdwtyp == "LENNARD-JONES") elj;
    //     if (vdwtyp == "BUCKINGHAM") ebuck;
    //     if (vdwtyp == "MM3-HBOND") emm3hb;
    //     if (vdwtyp == "BUFFERED-14-7") ehal;
    //     if (vdwtyp == "GAUSSIAN") egauss;
    // }
    // if (use_repel) erepel();
    // if (use_disp) edisp();

    // call any miscellaneous energy component routines
    // if (use_solv) esolv();
    // if (use_geom) egeom();
    // if (use_metal) emetal();
    // if (use_extra) extra();

    // sum up to give the total potential energy
    if constexpr (do_e) {
        esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
            + et + ept + ebt + eat + ett + ev + er + edsp
            + ec + ecd + ed + em + ep + ect + erxf + es + elf
            + eg + ex;
    }

    // sum up to give the total potential energy per atom
    if constexpr (do_a) {
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            aesum[i] = aeb[i] + aea[i] + aeba[i] + aeub[i] + aeaa[i]
                + aeopb[i] + aeopd[i] + aeid[i] + aeit[i]
                + aet[i] + aept[i] + aebt[i] + aeat[i] + aett[i]
                + aev[i] + aer[i] + aedsp[i] + aec[i] + aecd[i]
                + aed[i] + aem[i] + aep[i] + aect[i] + aerxf[i]
                + aes[i] + aelf[i] + aeg[i] + aex[i];
        }
    }

    // sum up to give the total derivative per atom
    if constexpr (do_g) {
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                desum[3*i+j] = deb[3*i+j] + dea[3*i+j] + deba[3*i+j]
                    + deub[3*i+j] + deaa[3*i+j] + deopb[3*i+j]
                    + deopd[3*i+j] + deid[3*i+j] + deit[3*i+j]
                    + det[3*i+j] + dept[3*i+j] + debt[3*i+j]
                    + deat[3*i+j] + dett[3*i+j] + dev[3*i+j]
                    + der[3*i+j] + dedsp[3*i+j] + dec[3*i+j]
                    + decd[3*i+j] + ded[3*i+j] + dem[3*i+j]
                    + dep[3*i+j] + dect[3*i+j] + derxf[3*i+j]
                    + des[3*i+j] + delf[3*i+j]
                    + deg[3*i+j] + dex[3*i+j];
            }
        }
    }

    // check for an illegal value for the total energy
    if (std::isnan(esum)) {
        printf("\n ENERGY  --  Illegal Value for the Total Potential Energy\n");
        fatal();
    }
}

// explicit instatiation
template void energy<CalcMode::Energy>();
template void energy<CalcMode::Analysis>();
template void energy<CalcMode::Gradient>();
template void energy<CalcMode::Virial>();

void initEnergy()
{
    // initial dynamic allocation of some global arrays
    aesum.allocate(1);
    aeb.allocate(1);
    aea.allocate(1);
    aeba.allocate(1);
    aeub.allocate(1);
    aeaa.allocate(1);
    aeopb.allocate(1);
    aeopd.allocate(1);
    aeid.allocate(1);
    aeit.allocate(1);
    aet.allocate(1);
    aept.allocate(1);
    aebt.allocate(1);
    aeat.allocate(1);
    aett.allocate(1);
    aev.allocate(1);
    aer.allocate(1);
    aedsp.allocate(1);
    aec.allocate(1);
    aecd.allocate(1);
    aed.allocate(1);
    aem.allocate(1);
    aep.allocate(1);
    aect.allocate(1);
    aerxf.allocate(1);
    aes.allocate(1);
    aelf.allocate(1);
    aeg.allocate(1);
    aex.allocate(1);
    desum.allocate(1);
    deb.allocate(1);
    dea.allocate(1);
    deba.allocate(1);
    deub.allocate(1);
    deaa.allocate(1);
    deopb.allocate(1);
    deopd.allocate(1);
    deid.allocate(1);
    deit.allocate(1);
    det.allocate(1);
    dept.allocate(1);
    debt.allocate(1);
    deat.allocate(1);
    dett.allocate(1);
    dev.allocate(1);
    der.allocate(1);
    dedsp.allocate(1);
    dec.allocate(1);
    decd.allocate(1);
    ded.allocate(1);
    dem.allocate(1);
    dep.allocate(1);
    dect.allocate(1);
    derxf.allocate(1);
    des.allocate(1);
    delf.allocate(1);
    deg.allocate(1);
    dex.allocate(1);
    te.allocate(1);
}
}
