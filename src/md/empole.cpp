// Author: Moses KJ Chung
// Year:   2023

#include "action.h"
#include "analyz.h"
#include "atomid.h"
#include "atoms.h"
#include "bound.h"
#include "calcMode.h"
#include "cell.h"
#include "chgpen.h"
#include "chgpot.h"
#include "chkpole.h"
#include "couple.h"
#include "cutoffSwitch.h"
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
            if (pentype == PenTyp::None)
                empole_a<CalculationMode, PenTyp::None>();
            else if (pentype == PenTyp::Gordon1)
                empole_a<CalculationMode, PenTyp::Gordon1>();
            else if (pentype == PenTyp::Gordon2)
                empole_a<CalculationMode, PenTyp::Gordon2>();
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

template <CalcMode CalculationMode, PenTyp PenType>
void empole_a()
{
    int i,k;
    int ix,iy,iz;
    int kx,ky,kz;
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
    real ttmxi,ttmyi,ttmzi;
    real ttmxk,ttmyk,ttmzk;
    real vxx,vyy,vzz;
    real vxy,vxz,vyz;
    std::vector<real> mscale;
    bool proceed;
    bool usei,usek;
    MAYBE_UNUSED std::vector<std::vector<real>> tem;

    if (npole == 0) return;

    // choose calculation mode
    constexpr CalcFlag flags = getCalculationFlags<CalculationMode>();
    constexpr bool do_e = flags.do_energy;
    constexpr bool do_a = flags.do_analysis;
    constexpr bool do_g = flags.do_gradient;
    constexpr bool do_v = flags.do_virial;

    // zero out total atomic multipole energy and partitioning
    em = 0.;
    if constexpr (do_a) {
        nem = 0;
        for (int i = 0; i < n; i++) {
            aem[i] = 0.;
        }
    }

    // check the sign of multipole components at chiral sites
    chkpole();

    // rotate the multipole components into the global frame
    rotpole(RotMode::Mpole);

    // allocate and initialize connected atom exclusion coefficients
    mscale.resize(n, 1.);

    // allocate and initialize torque array
    if constexpr (do_g) tem.resize(n, std::vector<real>(3,0.));

    // set conversion factor, cutoff and switching coefficients
    f = electric / dielec;
    cutoffSwitch(CutoffMode::Mpole);

    // calculate the multipole interaction energy term
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
                        pairMpoleCP_a<do_e, do_g, do_v, PenType>(
                            r2, xr, yr, zr, mk,
                            corei, vali, alphai, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz,
                            corek, valk, alphak, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,
                            f, frcx, frcy, frcz,
                            ttmxi, ttmyi, ttmzi, ttmxk, ttmyk, ttmzk,
                            e, vxx, vxy, vxz, vyy, vyz, vzz);
                    }
                    else {
                        pairMpole_a<do_e, do_g, do_v>(
                            r2, xr, yr, zr, mk,
                            ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz,
                            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,
                            f, frcx, frcy, frcz,
                            ttmxi, ttmyi, ttmzi, ttmxk, ttmyk, ttmzk,
                            e, vxx, vxy, vxz, vyy, vyz, vzz);
                    }

                    // increment the overall multipoles
                    if constexpr (do_a) {
                        if (e != 0.) nem++;
                        aem[i] = aem[i] + 0.5*e;
                        aem[k] = aem[k] + 0.5*e;
                        if (molcule[i] != molcule[k]) {
                            einter += e;
                        }
                    }
                    if constexpr (do_e) {
                        em += e;
                    }
                    if constexpr (do_g) {
                        dem[i][0] += frcx;
                        dem[i][1] += frcy;
                        dem[i][2] += frcz;
                        tem[i][0] += ttmxi;
                        tem[i][1] += ttmyi;
                        tem[i][2] += ttmzi;
                        dem[k][0] -= frcx;
                        dem[k][1] -= frcy;
                        dem[k][2] -= frcz;
                        tem[k][0] += ttmxk;
                        tem[k][1] += ttmyk;
                        tem[k][2] += ttmzk;
                    }
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

    // resolve site torques then increment forces and virial
    if constexpr (do_g or do_v) {
        torque<CalculationMode>(&tem, &dem);
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
