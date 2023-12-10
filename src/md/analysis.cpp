// Author: Moses KJ Chung
// Year:   2023

#include "analysis.h"
#include "analyz.h"
#include "atoms.h"
#include "calcMode.h"
#include "empole.h"
#include "energi.h"
#include "fatal.h"
#include "group.h"
#include "inter.h"
#include "limits.h"
#include "potent.h"
#include "vdwpot.h"
#include <cmath>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  analysis  --  energy components and analysis  //
//                                                //
////////////////////////////////////////////////////

// "analysis" calls the series of routines needed to calculate
// the potential energy and perform energy partitioning analysis
// in terms of type of interaction or atom number

void analysis(double& energy)
{
    int i;
    double cutoff;

    calcMode calculationMode = calcMode::ANALYSIS;

    // zero out each of the potential energy components
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

    // perform dynamic allocation of some global arrays
    if (aesum.size() != 0) {
        if (aesum.size() < n) {
            aesum.resize(0);
            aeb.resize(0);
            aea.resize(0);
            aeba.resize(0);
            aeub.resize(0);
            aeaa.resize(0);
            aeopb.resize(0);
            aeopd.resize(0);
            aeid.resize(0);
            aeit.resize(0);
            aet.resize(0);
            aept.resize(0);
            aebt.resize(0);
            aeat.resize(0);
            aett.resize(0);
            aev.resize(0);
            aer.resize(0);
            aedsp.resize(0);
            aec.resize(0);
            aecd.resize(0);
            aed.resize(0);
            aem.resize(0);
            aep.resize(0);
            aect.resize(0);
            aerxf.resize(0);
            aes.resize(0);
            aelf.resize(0);
            aeg.resize(0);
            aex.resize(0);
        }
    }

    if (aesum.size() == 0) {
        aesum.resize(n);
        aeb.resize(n);
        aea.resize(n);
        aeba.resize(n);
        aeub.resize(n);
        aeaa.resize(n);
        aeopb.resize(n);
        aeopd.resize(n);
        aeid.resize(n);
        aeit.resize(n);
        aet.resize(n);
        aept.resize(n);
        aebt.resize(n);
        aeat.resize(n);
        aett.resize(n);
        aev.resize(n);
        aer.resize(n);
        aedsp.resize(n);
        aec.resize(n);
        aecd.resize(n);
        aed.resize(n);
        aem.resize(n);
        aep.resize(n);
        aect.resize(n);
        aerxf.resize(n);
        aes.resize(n);
        aelf.resize(n);
        aeg.resize(n);
        aex.resize(n);
    }

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

    // zero out the total intermolecular energy
    einter = 0.;

    // remove any previous use of the replicates method
    // cutoff = 0.;
    // replica(cutoff);

    // update the pairwise interaction neighbor lists
    // if (use_list) nblist(); // TODO

    // many implicit solvation models require Born radii
    // if (use_born) born(); // TODO

    // alter partial charges and multipoles for charge flux
    // if (use_chgflx) alterchg(); // TODO

    // modify bond and torsion constants for pisystem
    // if (use_orbit) picalc(); // TODO

    // call the local geometry energy component routines
    // if (use_bond) ebond3(); // TODO
    // if (use_angle) eangle3(); // TODO
    // if (use_strbnd) estrbnd3(); // TODO
    // if (use_urey) eurey3(); // TODO
    // if (use_angang) eangang3(); // TODO
    // if (use_opbend) eopbend3(); // TODO
    // if (use_opdist) eopdist3(); // TODO
    // if (use_improp) eimprop3(); // TODO
    // if (use_imptor) eimptor3(); // TODO
    // if (use_tors) etors3(); // TODO
    // if (use_pitors) epitors3(); // TODO
    // if (use_strtor) estrtor3(); // TODO
    // if (use_angtor) eangtor3(); // TODO
    // if (use_tortor) etortor3(); // TODO

    // call the electrostatic energy component routines
    // if (use_charge) echarge3(); // TODO
    // if (use_chgdpl) echgdpl3(); // TODO
    // if (use_dipole) edipole3(); // TODO
    if (use_mpole) empole(calculationMode);
    // if (use_polar) epolar3(); // TODO
    // if (use_chgtrn) echgtrn3(); // TODO
    // if (use_rxnfld) erxnfld3(); // TODO

    // call the van der Waals energy component routines
    if (use_vdw) {
        // if (vdwtyp == "LENNARD-JONES") elj3(); // TODO
        // if (vdwtyp == "BUCKINGHAM") ebuck3(); // TODO
        // if (vdwtyp == "MM3-HBOND") emm3hb3(); // TODO
        // if (vdwtyp == "BUFFERED-14-7") ehal3(); // TODO
        // if (vdwtyp == "GAUSSIAN") egauss3(); // TODO
    }
    // if (use_repel) erepel3(); // TODO
    // if (use_disp) edisp3(); // TODO

    // call any miscellaneous energy component routines
    // if (use_solv) esolv3(); // TODO
    // if (use_metal) emetal3(); // TODO
    // if (use_geom) egeom3(); // TODO
    // if (use_extra) extra3(); // TODO

    // sum up to give the total potential energy
    esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
              + et + ept + ebt + eat + ett + ev + er + edsp
              + ec + ecd + ed + em + ep + ect + erxf + es + elf
              + eg + ex;
    energy = esum;

    // sum up to give the total potential energy per atom
    for (int i = 0; i < n; i++) {
        aesum[i] = aeb[i] + aea[i] + aeba[i] + aeub[i] + aeaa[i]
                          + aeopb[i] + aeopd[i] + aeid[i] + aeit[i]
                          + aet[i] + aept[i] + aebt[i] + aeat[i] + aett[i]
                          + aev[i] + aer[i] + aedsp[i] + aec[i] + aecd[i]
                          + aed[i] + aem[i] + aep[i] + aect[i] + aerxf[i]
                          + aes[i] + aelf[i] + aeg[i] + aex[i];
    }

    // check for an illegal value for the total energy
    if (std::isnan(esum)) {
        printf("\n ANALYSIS  --  Illegal Value for the Total Potential Energy\n");
        fatal();
    }
}
}
