// Author: Moses KJ Chung
// Year:   2024

#include "action.h"
#include "analyz.h"
#include "atoms.h"
#include "bound.h"
#include "chkpole.h"
#include "couple.h"
#include "cutoffSwitch.h"
#include "deriv.h"
#include "energi.h"
#include "erepel.h"
#include "inter.h"
#include "mdqclimits.h"
#include "molcul.h"
#include "potent.h"
#include "repel.h"
#include "reppot.h"
#include "rotpole.h"
#include "shunt.h"
#include "usage.h"
#include "virial.h"

namespace polmdqc
{
///////////////////////////////////////////////
//                                           //
//  erepel  --  Pauli repulsion calculation  //
//                                           //
///////////////////////////////////////////////

// "erepel" calculates the Pauli repulsion energy
// and/or gradient among the atoms
//
// literature reference:
//
// J. A. Rackers and J. W. Ponder, "Classical Pauli Repulsion:
// An Anisotropic, Atomic Multipole Model", Journal of Chemical
// Physics, 150, 084104 (2019)

template <CalcMode CalculationMode>
void erepel()
{
    // choose the method for summing over pairwise interactions
    if (use_mlist) {
        // erepel_b<CalculationMode>();
    }
    else {
        erepel_a<CalculationMode>();
    }
}

/////////////////////////////////////////////////////////////
//                                                         //
//  erepel_a  --  double loop Pauli repulsion calculation  //
//                                                         //
/////////////////////////////////////////////////////////////

// "erepel_a" calculates the Pauli repulsion interactions
// using a double loop

template <CalcMode CalculationMode>
void erepel_a()
{
    int tid;
    real mk;
    real e;
    real xi,yi,zi;
    real xr,yr,zr;
    real r2;
    real ci,dix,diy,diz;
    real qixx,qixy,qixz;
    real qiyy,qiyz,qizz;
    real ck,dkx,dky,dkz;
    real qkxx,qkxy,qkxz;
    real qkyy,qkyz,qkzz;
    real sizi,sizk,sizik;
    real vali,valk;
    real dmpi,dmpk;
    real frcx,frcy,frcz;
    real ttrxi,ttryi,ttrzi;
    real ttrxk,ttryk,ttrzk;
    real vxx,vyy,vzz;
    real vxy,vxz,vyz;
    real* rscale;

    if (nrep == 0) return;

    // choose calculation mode
    constexpr CalcFlag flags = getCalculationFlags<CalculationMode>();
    constexpr bool do_e = flags.do_energy;
    constexpr bool do_a = flags.do_analysis;
    constexpr bool do_g = flags.do_gradient;
    constexpr bool do_v = flags.do_virial;

    // set pointers for OpenMP
    real* aerP = aer.ptr();
    real* derP = der.ptr();
    real* terP = te.ptr();

    // zero out total repulsion energy and partitioning
    er = 0.;
    if constexpr (do_a) {
        ner = 0;
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            aerP[i] = 0.;
        }
    }

    // zero out torque array
    if constexpr (do_g) {
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                terP[3*i+j] = 0.;
            }
        }
    }

    // check the sign of multipole components at chiral sites
    chkrepole();

    // rotate the multipole components into the global frame
    rotpole(RotMode::Repel);

    // set conversion factor, cutoff and switching coefficients
    cutoffSwitch(CutoffMode::Repuls);

    // OpenMP setup
    int Ndo_a = 1;
    int Ndo_3g = 1;
    if (do_a) Ndo_a = n;
    if (do_g) Ndo_3g = 3*n;
    #pragma omp parallel default(private)                                                 \
    shared(n,xaxis,yaxis,zaxis,x,y,z,rrepole,sizpr,dmppr,elepr,use,use_bounds,off2,       \
        n12,i12,n13,i13,n14,i14,n15,i15,r2scale,r3scale,r4scale,r5scale,molcule,escale)   \
    reduction(+:er,einter,ner,vir,aerP[:Ndo_a],derP[:Ndo_3g],terP[:Ndo_3g])
    {
    // initialize connected atom exclusion coefficients
    tid = omp_get_thread_num();
    rscale = escale[tid];
    for (int i = 0; i < n; i++) {
        rscale[i] = 1.;
    }

    // calculate the Pauli repulsion interaction term
    #pragma omp for schedule(guided)
    for (int i = 0; i < n-1; i++) {
    }
    }
}

// explicit instatiation
template void erepel<CalcMode::Energy>();
template void erepel<CalcMode::Analysis>();
template void erepel<CalcMode::Gradient>();
template void erepel<CalcMode::Virial>();
}
