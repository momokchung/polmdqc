// Author: Moses KJ Chung
// Year:   2023

#include "cutoffSwitch.h"
#include "nonpol.h"
#include "mdqclimits.h"
#include "shunt.h"
#include <algorithm>
#include <cmath>

namespace polmdqc
{
/////////////////////////////////////////////////////////////
//                                                         //
//  cutoffSwitch  --  get switching function coefficients  //
//                                                         //
/////////////////////////////////////////////////////////////

// "cutoffSwitch" sets the coeffcients used by the fifth and seventh
// order polynomial switching functions for spherical cutoffs

void cutoffSwitch(CutoffMode mode)
{
    double denom,term;
    double off3,off4,off5;
    double off6,off7;
    double cut3,cut4,cut5;
    double cut6,cut7;

    // get the switching window for the current potential type
    if (mode == CutoffMode::VdW) {
        off = vdwcut;
        cut = vdwtaper;
    }
    else if (mode == CutoffMode::Repuls) {
        off = repcut;
        cut = reptaper;
    }
    else if (mode == CutoffMode::Disp) {
        off = dispcut;
        cut = disptaper;
    }
    else if (mode == CutoffMode::Charge) {
        off = chgcut;
        cut = chgtaper;
    }
    else if (mode == CutoffMode::ChgDpl) {
        off = std::sqrt(chgcut*dplcut);
        cut = std::sqrt(chgtaper*dpltaper);
    }
    else if (mode == CutoffMode::Dipole) {
        off = dplcut;
        cut = dpltaper;
    }
    else if (mode == CutoffMode::Mpole) {
        off = mpolecut;
        cut = mpoletaper;
    }
    else if (mode == CutoffMode::ChgTrn) {
        off = ctrncut;
        cut = ctrntaper;
    }
    else if (mode == CutoffMode::Ewald) {
        off = ewaldcut;
        cut = ewaldcut;
    }
    else if (mode == CutoffMode::DEwald) {
        off = dewaldcut;
        cut = dewaldcut;
    }
    else if (mode == CutoffMode::USolv) {
        off = usolvcut;
        cut = usolvcut;
    }
    else if (mode == CutoffMode::GKV) {
        off = spoff;
        cut = spcut;
    }
    else if (mode == CutoffMode::GKSA) {
        off = stcut;
        cut = stoff;
    }
    else {
        off = std::min({vdwcut,repcut,dispcut,chgcut,dplcut,mpolecut,ctrncut});
        cut = std::min({vdwtaper,reptaper,disptaper,chgtaper,dpltaper,mpoletaper,ctrntaper});
    }

    // test for replicate periodic boundaries at this cutoff
    // replica(off);

    // set switching coefficients to zero for truncation cutoffs
    c0 = 0.;
    c1 = 0.;
    c2 = 0.;
    c3 = 0.;
    c4 = 0.;
    c5 = 0.;
    f0 = 0.;
    f1 = 0.;
    f2 = 0.;
    f3 = 0.;
    f4 = 0.;
    f5 = 0.;
    f6 = 0.;
    f7 = 0.;

    // store the powers of the switching window cutoffs
    off2 = off * off;
    off3 = off2 * off;
    off4 = off2 * off2;
    off5 = off2 * off3;
    off6 = off3 * off3;
    off7 = off3 * off4;
    cut2 = cut * cut;
    cut3 = cut2 * cut;
    cut4 = cut2 * cut2;
    cut5 = cut2 * cut3;
    cut6 = cut3 * cut3;
    cut7 = cut3 * cut4;

    // get 5th degree multiplicative switching function coefficients
    if (cut < off) {
        denom = std::pow((off-cut),5);
        c0 = off*off2 * (off2-5.*off*cut+10.*cut2) / denom;
        c1 = -30. * off2*cut2 / denom;
        c2 = 30. * (off2*cut+off*cut2) / denom;
        c3 = -10. * (off2+4.*off*cut+cut2) / denom;
        c4 = 15. * (off+cut) / denom;
        c5 = -6. / denom;
    }

    // get 7th degree additive switching function coefficients
    if (cut<off and mode==CutoffMode::Charge) {
        term = 9.3 * cut*off / (off-cut);
        denom = cut7 - 7.*cut6*off + 21.*cut5*off2 - 35.*cut4*off3 + 35.*cut3*off4 - 21.*cut2*off5 + 7.*cut*off6 - off7;
        denom = term * denom;
        f0 = cut3*off3 * (-39.*cut+64.*off) / denom;
        f1 = cut2*off2 * (117.*cut2-100.*cut*off-192.*off2) / denom;
        f2 = cut*off * (-117.*cut3-84.*cut2*off+534.*cut*off2+192.*off3) / denom;
        f3 = (39.*cut4+212.*cut3*off-450.*cut2*off2-612.*cut*off3-64.*off4) / denom;
        f4 = (-92.*cut3+66.*cut2*off+684.*cut*off2+217.*off3) / denom;
        f5 = (42.*cut2-300.*cut*off-267.*off2) / denom;
        f6 = (36.*cut+139.*off) / denom;
        f7 = -25. / denom;
    }
}
}
