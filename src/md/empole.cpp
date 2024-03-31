// Author: Moses KJ Chung
// Year:   2023

#include "action.h"
#include "analyz.h"
#include "atomid.h"
#include "atoms.h"
#include "bound.h"
#include "calcMode.h"
#include "cell.h"
#include "cflux.h"
#include "chgpen.h"
#include "chgpot.h"
#include "chkpole.h"
#include "couple.h"
#include "cutoffSwitch.h"
#include "dcflux.h"
#include "deriv.h"
#include "empole.h"
#include "energi.h"
#include "image.h"
#include "inform.h"
#include "inter.h"
#include "mathConst.h"
#include "molcul.h"
#include "mpole.h"
#include "pairMpole.h"
#include "potent.h"
#include "mdqclimits.h"
#include "rotpole.h"
#include "shunt.h"
#include "torque.h"
#include "usage.h"
#include "virial.h"
#include <cmath>

namespace polmdqc
{
////////////////////////////////////////////////
//                                            //
//  empole  --  atomic multipole calculation  //
//                                            //
////////////////////////////////////////////////

// "empole" calculates the electrostatic energy and/or
// gradient due to atomic multipole interactions

template <CalcMode CalculationMode>
void empole()
{
    // choose the method to sum over multipole interactions
    if (use_ewald) {
        if (use_mlist) {
            // empole_d<CalculationMode>();
        }
        else {
            // empole_c<CalculationMode>();
        }
    }
    else {
        if (use_mlist) {
            // empole_b<CalculationMode>();
        }
        else {
            if (pentype == PenTyp::None) {
                empole_a<CalculationMode, PenTyp::None, false>();
            }
            else if (pentype == PenTyp::Gordon1) {
                if (use_chgflx)
                    empole_a<CalculationMode, PenTyp::Gordon1, true>();
                else
                    empole_a<CalculationMode, PenTyp::Gordon1, false>();
            }
            else if (pentype == PenTyp::Gordon2) {
                if (use_chgflx)
                    empole_a<CalculationMode, PenTyp::Gordon2, true>();
                else
                    empole_a<CalculationMode, PenTyp::Gordon2, false>();
            }
        }
    }
}

///////////////////////////////////////////////////////
//                                                   //
//  empole_a  --  double loop multipole calculation  //
//                                                   //
///////////////////////////////////////////////////////

// "empole_a" calculates the atomic multipole interactions
// using a double loop

template <CalcMode CalculationMode, PenTyp PenType, bool use_cf>
void empole_a()
{
    int ix,iy,iz;
    int kx,ky,kz;
    int tid;
    real mk;
    real e,f;
    real xi,yi,zi;
    real xr,yr,zr;
    real r2;
    real corei,vali,alphai;
    real ci,dix,diy,diz;
    real qixx,qixy,qixz;
    real qiyy,qiyz,qizz;
    real corek,valk,alphak;
    real ck,dkx,dky,dkz;
    real qkxx,qkxy,qkxz;
    real qkyy,qkyz,qkzz;
    real frcx,frcy,frcz;
    real trqxi,trqyi,trqzi;
    real trqxk,trqyk,trqzk;
    real vxx,vyy,vzz;
    real vxy,vxz,vyz;
    real poti,potk;
    real* mscale;
    bool proceed;
    bool usei,usek;

    if (npole == 0) return;

    // choose calculation mode
    constexpr CalcFlag flags = getCalculationFlags<CalculationMode>();
    constexpr bool do_e = flags.do_energy;
    constexpr bool do_a = flags.do_analysis;
    constexpr bool do_g = flags.do_gradient;
    constexpr bool do_v = flags.do_virial;

    // set pointers for OpenMP
    real* aemP = aem.ptr();
    real* demP = dem.ptr();
    real* temP = te.ptr();
    real* potP = pot.ptr();

    // zero out total atomic multipole energy and partitioning
    em = 0.;
    if constexpr (do_a) {
        nem = 0;
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            aemP[i] = 0.;
        }
    }

    // zero out charge flux potential array
    if constexpr (do_g and use_cf) {
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            potP[i] = 0.;
        }
    }

    // zero out torque array
    if constexpr (do_g) {
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                temP[3*i+j] = 0.;
            }
        }
    }

    // check the sign of multipole components at chiral sites
    chkpole();

    // rotate the multipole components into the global frame
    rotpole(RotMode::Mpole);

    // set conversion factor, cutoff and switching coefficients
    f = electric / dielec;
    cutoffSwitch(CutoffMode::Mpole);

    // OpenMP setup
    int Ndo_a = 1;
    int Ndo_g = 1;
    int Ndo_3g = 1;
    if (do_a) Ndo_a = n;
    if (do_g and use_cf) Ndo_g = n;
    if (do_g) Ndo_3g = 3*n;
    #pragma omp parallel default(private)                                                 \
    shared(n,xaxis,yaxis,zaxis,x,y,z,rpole,pcore,pval,palpha,use,use_bounds,off2,f,       \
        n12,i12,n13,i13,n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,molcule,escale)   \
    reduction(+:em,einter,nem,vir,aemP[:Ndo_a],potP[:n],demP[:Ndo_3g],temP[:Ndo_3g])
    {
    // initialize connected atom exclusion coefficients
    tid = omp_get_thread_num();
    mscale = escale[tid];
    for (int i = 0; i < n; i++) {
        mscale[i] = 1.;
    }

    // calculate the multipole interaction term
    #pragma omp for schedule(guided)
    for (int i = 0; i < n-1; i++) {
        iz = zaxis[i] - 1;
        ix = xaxis[i] - 1;
        iy = std::abs(yaxis[i]) - 1;
        xi = x[i];
        yi = y[i];
        zi = z[i];
        ci = rpole[i][0];
        dix = rpole[i][1];
        diy = rpole[i][2];
        diz = rpole[i][3];
        qixx = rpole[i][4];
        qixy = rpole[i][5];
        qixz = rpole[i][6];
        qiyy = rpole[i][8];
        qiyz = rpole[i][9];
        qizz = rpole[i][12];
        if constexpr (PenType == PenTyp::Gordon1 or PenType == PenTyp::Gordon2) {
            corei = pcore[i];
            vali = pval[i];
            alphai = palpha[i];
        }
        usei = (use[i+1] or use[iz+1] or use[ix+1] or use[iy+1]);

        // set exclusion coefficients for connected atoms
        for (int j = 0; j < n12[i]; j++) {
            mscale[i12[i][j]] = m2scale;
        }
        for (int j = 0; j < n13[i]; j++) {
            mscale[i13[i][j]] = m3scale;
        }
        for (int j = 0; j < n14[i]; j++) {
            mscale[i14[i][j]] = m4scale;
        }
        for (int j = 0; j < n15[i]; j++) {
            mscale[i15[i][j]] = m5scale;
        }

        // evaluate all sites within the cutoff distance
        for (int k = i+1; k < n; k++) {
            kz = zaxis[k] - 1;
            kx = xaxis[k] - 1;
            ky = std::abs(yaxis[k]) - 1;
            usek = (use[k+1] or use[kz+1] or use[kx+1] or use[ky+1]);
            proceed = (usei or usek);
            if (proceed) {
                xr = x[k] - xi;
                yr = y[k] - yi;
                zr = z[k] - zi;
                if (use_bounds) image(xr,yr,zr);
                r2 = xr * xr + yr * yr + zr * zr;
                if (r2 <= off2) {
                    mk = mscale[k];
                    ck = rpole[k][0];
                    dkx = rpole[k][1];
                    dky = rpole[k][2];
                    dkz = rpole[k][3];
                    qkxx = rpole[k][4];
                    qkxy = rpole[k][5];
                    qkxz = rpole[k][6];
                    qkyy = rpole[k][8];
                    qkyz = rpole[k][9];
                    qkzz = rpole[k][12];
                    if constexpr (PenType == PenTyp::Gordon1 or PenType == PenTyp::Gordon2) {
                        corek = pcore[k];
                        valk = pval[k];
                        alphak = palpha[k];
                        pairMpoleCP<do_e, do_g, do_v, PenType, use_cf>(
                            r2, xr, yr, zr, mk,
                            corei, vali, alphai, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz,
                            corek, valk, alphak, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,
                            f, frcx, frcy, frcz,
                            trqxi, trqyi, trqzi, trqxk, trqyk, trqzk,
                            e, vxx, vxy, vxz, vyy, vyz, vzz, poti, potk);
                    }
                    else {
                        pairMpole<do_e, do_g, do_v>(
                            r2, xr, yr, zr, mk,
                            ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz,
                            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,
                            f, frcx, frcy, frcz,
                            trqxi, trqyi, trqzi, trqxk, trqyk, trqzk,
                            e, vxx, vxy, vxz, vyy, vyz, vzz);
                    }

                    // increment the overall multipoles
                    if constexpr (do_a) {
                        if (e != 0.) nem++;
                        aemP[i] += 0.5*e;
                        aemP[k] += 0.5*e;
                        if (molcule[i] != molcule[k]) {
                            einter += e;
                        }
                    }
                    if constexpr (do_e) {
                        em += e;
                    }
                    if constexpr (do_g) {
                        demP[3*i+0] += frcx;
                        demP[3*i+1] += frcy;
                        demP[3*i+2] += frcz;
                        temP[3*i+0] += trqxi;
                        temP[3*i+1] += trqyi;
                        temP[3*i+2] += trqzi;
                        demP[3*k+0] -= frcx;
                        demP[3*k+1] -= frcy;
                        demP[3*k+2] -= frcz;
                        temP[3*k+0] += trqxk;
                        temP[3*k+1] += trqyk;
                        temP[3*k+2] += trqzk;

                        if constexpr (do_v) {
                            vir[0][0] += vxx;
                            vir[0][1] += vxy;
                            vir[0][2] += vxz;
                            vir[1][0] += vxy;
                            vir[1][1] += vyy;
                            vir[1][2] += vyz;
                            vir[2][0] += vxz;
                            vir[2][1] += vyz;
                            vir[2][2] += vzz;
                        }

                        if constexpr (use_cf) {
                            potP[i] += poti;
                            potP[k] += potk;
                        }
                    }
                }
            }
        }

        // reset exclusion coefficients for connected atoms
        for (int j = 0; j < n12[i]; j++) {
            mscale[i12[i][j]] = 1.;
        }
        for (int j = 0; j < n13[i]; j++) {
            mscale[i13[i][j]] = 1.;
        }
        for (int j = 0; j < n14[i]; j++) {
            mscale[i14[i][j]] = 1.;
        }
        for (int j = 0; j < n15[i]; j++) {
            mscale[i15[i][j]] = 1.;
        }

    }
    }

    // resolve site torques then increment forces and virial
    if constexpr (do_g) {
        torque<CalculationMode>(temP, demP, vir);
    }

    // modify the gradient and virial for charge flux
    if constexpr (do_g and use_cf) {
        dcflux<CalculationMode>(potP, demP, vir);
    }
}

// for neighborlist, use Zhi's method of dividing calculations
// into exclusion/nonexclusion pairs
// template <CalcMode CalculationMode>
// void empole_b(){}

// template <CalcMode CalculationMode>
// void empole_c(){}

// template <CalcMode CalculationMode>
// void empole_d(){}

// explicit instatiation
template void empole<CalcMode::Energy>();
template void empole<CalcMode::Analysis>();
template void empole<CalcMode::Gradient>();
template void empole<CalcMode::Virial>();
}
