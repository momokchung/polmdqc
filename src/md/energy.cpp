// Author: Moses KJ Chung
// Year:   2024

#include "analyz.h"
#include "atoms.h"
#include "calcMode.h"
#include "deriv.h"
#include "empole.h"
#include "energi.h"
#include "energy.h"
#include "fatal.h"
#include "inter.h"
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

        // zero out the total intermolecular energy
        einter = 0.;
    }

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
        desum.allocate(n);
        deb.allocate(n);
        dea.allocate(n);
        deba.allocate(n);
        deub.allocate(n);
        deaa.allocate(n);
        deopb.allocate(n);
        deopd.allocate(n);
        deid.allocate(n);
        deit.allocate(n);
        det.allocate(n);
        dept.allocate(n);
        debt.allocate(n);
        deat.allocate(n);
        dett.allocate(n);
        dev.allocate(n);
        der.allocate(n);
        dedsp.allocate(n);
        dec.allocate(n);
        decd.allocate(n);
        ded.allocate(n);
        dem.allocate(n);
        dep.allocate(n);
        dect.allocate(n);
        derxf.allocate(n);
        des.allocate(n);
        delf.allocate(n);
        deg.allocate(n);
        dex.allocate(n);

        // zero out each of the first derivative components
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                desum[i][j] = 0.;
                deb[i][j] = 0.;
                dea[i][j] = 0.;
                deba[i][j] = 0.;
                deub[i][j] = 0.;
                deaa[i][j] = 0.;
                deopb[i][j] = 0.;
                deopd[i][j] = 0.;
                deid[i][j] = 0.;
                deit[i][j] = 0.;
                det[i][j] = 0.;
                dept[i][j] = 0.;
                debt[i][j] = 0.;
                deat[i][j] = 0.;
                dett[i][j] = 0.;
                dev[i][j] = 0.;
                der[i][j] = 0.;
                dedsp[i][j] = 0.;
                dec[i][j] = 0.;
                decd[i][j] = 0.;
                ded[i][j] = 0.;
                dem[i][j] = 0.;
                dep[i][j] = 0.;
                dect[i][j] = 0.;
                derxf[i][j] = 0.;
                des[i][j] = 0.;
                delf[i][j] = 0.;
                deg[i][j] = 0.;
                dex[i][j] = 0.;
            }
        }
    }

    if constexpr (do_v) {
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
    // if (use_chgflx) alterchg();

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
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                desum[i][j] = deb[i][j] + dea[i][j] + deba[i][j]
                    + deub[i][j] + deaa[i][j] + deopb[i][j]
                    + deopd[i][j] + deid[i][j] + deit[i][j]
                    + det[i][j] + dept[i][j] + debt[i][j]
                    + deat[i][j] + dett[i][j] + dev[i][j]
                    + der[i][j] + dedsp[i][j] + dec[i][j]
                    + decd[i][j] + ded[i][j] + dem[i][j]
                    + dep[i][j] + dect[i][j] + derxf[i][j]
                    + des[i][j] + delf[i][j]
                    + deg[i][j] + dex[i][j];
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
}
