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
#include "damping.h"
#include "empole.h"
#include "energi.h"
#include "group.h"
#include "groups.h"
#include "image.h"
#include "inform.h"
#include "inter.h"
#include "mathConst.h"
#include "molcul.h"
#include "mplpot.h"
#include "mpole.h"
#include "potent.h"
#include "mdqclimits.h"
#include "rotpole.h"
#include "shunt.h"
#include "usage.h"
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

void empole(calcMode calculationMode)
{
    // choose the method to sum over multipole interactions
    if (use_ewald) {
        if (use_mlist) {
            empole_d(calculationMode);
        }
        else {
            empole_c(calculationMode);
        }
    }
    else {
        if (use_mlist) {
            empole_b(calculationMode);
        }
        else {
            empole_a(calculationMode);
        }
    }
}

///////////////////////////////////////////////////////
//                                                   //
//  empole_a  --  double loop multipole calculation  //
//                                                   //
///////////////////////////////////////////////////////

// "empole3a" calculates the atomic multipole interaction energy
// and/or graident using a double loop

void empole_a(calcMode calculationMode)
{
    int i,j,k;
    int ii,kk;
    int ix,iy,iz;
    int kx,ky,kz;
    double e,f,fgrp;
    double xi,yi,zi;
    double xr,yr,zr;
    double r,r2,rr1,rr3;
    double rr5,rr7,rr9;
    double rr1i,rr3i,rr5i;
    double rr1k,rr3k,rr5k;
    double rr1ik,rr3ik,rr5ik;
    double rr7ik,rr9ik;
    double ci,dix,diy,diz;
    double qixx,qixy,qixz;
    double qiyy,qiyz,qizz;
    double ck,dkx,dky,dkz;
    double qkxx,qkxy,qkxz;
    double qkyy,qkyz,qkzz;
    double dir,dkr,dik,qik;
    double qix,qiy,qiz,qir;
    double qkx,qky,qkz,qkr;
    double diqk,dkqi,qiqk;
    double corei,corek;
    double vali,valk;
    double alphai,alphak;
    double term1,term2,term3;
    double term4,term5;
    double term1i,term2i,term3i;
    double term1k,term2k,term3k;
    double term1ik,term2ik,term3ik;
    double term4ik,term5ik;
    double dmpi[9],dmpk[9];
    double dmpik[9];
    std::vector<double> mscale;
    bool proceed;
    bool header,huge;
    bool usei,usek;
    std::string mode;
    bool do_a = false;
    bool do_g = false;
    if (calculationMode == calcMode::ANALYSIS) do_a = true;
    if (calculationMode == calcMode::GRADIENT) do_g = true;

    // zero out total atomic multipole energy and partitioning
    em = 0.;
    if (do_a) {
        nem = 0;
        for (int i = 0; i < n; i++) {
            aem[i] = 0.;
        }
    }
    if (npole == 0)  return

    // check the sign of multipole components at chiral sites
    chkpole();

    // rotate the multipole components into the global frame
    rotpole("MPOLE");

    // allocate and initialize connected atom exclusion coefficients
    mscale.resize(n, 1.);

    // set conversion factor, cutoff and switching coefficients
    f = electric / dielec;
    mode = "MPOLE";
    cutoffSwitch(mode);

    // print header information if debug output was requested
    header = true;
    if (debug and npole!=0) {
        header = false;
        printf("\n Individual Atomic Multipole Interactions :\n\n");
        printf(" Type              Atom Names               Distance        Energy\n\n");
    }

    // calculate the multipole interaction energy term
    for (int ii = 0; ii < npole-1; ii++) {
        i = ipole[ii];
        iz = zaxis[i] -1;
        ix = xaxis[i] -1;
        iy = std::abs(yaxis[i]) -1;
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
        if (use_chgpen) {
            corei = pcore[i];
            vali = pval[i];
            alphai = palpha[i];
        }
        usei = (use[i] or use[iz] or use[ix] or use[iy]);

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
        for (int kk = ii+1; kk < npole; kk++) {
            k = ipole[kk];
            kz = zaxis[k] - 1;
            kx = xaxis[k] - 1;
            ky = std::abs(yaxis[k]) - 1;
            usek = (use[k] or use[kz] or use[kx] or use[ky]);
            proceed = true;
            if (use_group) groups(proceed,fgrp,i,k,-1,-1,-1,-1);
            if (!use_intra) proceed = true;
            if (proceed) proceed = (usei or usek);
            if (proceed) {
                xr = x[k] - xi;
                yr = y[k] - yi;
                zr = z[k] - zi;
                if (use_bounds) image(xr,yr,zr);
                r2 = xr*xr + yr* yr + zr*zr;
                if (r2 <= off2) {
                    r = std::sqrt(r2);
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

                    // intermediates involving moments and separation distance
                    dir = dix*xr + diy*yr + diz*zr;
                    qix = qixx*xr + qixy*yr + qixz*zr;
                    qiy = qixy*xr + qiyy*yr + qiyz*zr;
                    qiz = qixz*xr + qiyz*yr + qizz*zr;
                    qir = qix*xr + qiy*yr + qiz*zr;
                    dkr = dkx*xr + dky*yr + dkz*zr;
                    qkx = qkxx*xr + qkxy*yr + qkxz*zr;
                    qky = qkxy*xr + qkyy*yr + qkyz*zr;
                    qkz = qkxz*xr + qkyz*yr + qkzz*zr;
                    qkr = qkx*xr + qky*yr + qkz*zr;
                    dik = dix*dkx + diy*dky + diz*dkz;
                    qik = qix*qkx + qiy*qky + qiz*qkz;
                    diqk = dix*qkx + diy*qky + diz*qkz;
                    dkqi = dkx*qix + dky*qiy + dkz*qiz;
                    qiqk = 2.*(qixy*qkxy+qixz*qkxz+qiyz*qkyz) + qixx*qkxx + qiyy*qkyy + qizz*qkzz;

                    // get reciprocal distance terms for this interaction
                    rr1 = f * mscale[k] / r;
                    rr3 = rr1 / r2;
                    rr5 = 3. * rr3 / r2;
                    rr7 = 5. * rr5 / r2;
                    rr9 = 7. * rr7 / r2;

                    // find damped multipole intermediates and energy value
                    if (use_chgpen) {
                        corek = pcore[k];
                        valk = pval[k];
                        alphak = palpha[k];
                        term1 = corei*corek;
                        term1i = corek*vali;
                        term2i = corek*dir;
                        term3i = corek*qir;
                        term1k = corei*valk;
                        term2k = -corei*dkr;
                        term3k = corei*qkr;
                        term1ik = vali*valk;
                        term2ik = valk*dir - vali*dkr + dik;
                        term3ik = vali*qkr + valk*qir - dir*dkr + 2.*(dkqi-diqk+qiqk);
                        term4ik = dir*qkr - dkr*qir - 4.*qik;
                        term5ik = qir*qkr;
                        damppole(r,9,alphai,alphak,dmpi,dmpk,dmpik);
                        rr1i = dmpi[0]*rr1;
                        rr3i = dmpi[2]*rr3;
                        rr5i = dmpi[4]*rr5;
                        rr1k = dmpk[0]*rr1;
                        rr3k = dmpk[2]*rr3;
                        rr5k = dmpk[4]*rr5;
                        rr1ik = dmpik[0]*rr1;
                        rr3ik = dmpik[2]*rr3;
                        rr5ik = dmpik[4]*rr5;
                        rr7ik = dmpik[6]*rr7;
                        rr9ik = dmpik[8]*rr9;
                        e = term1*rr1 + term1i*rr1i
                          + term1k*rr1k + term1ik*rr1ik
                          + term2i*rr3i + term2k*rr3k
                          + term2ik*rr3ik + term3i*rr5i
                          + term3k*rr5k + term3ik*rr5ik
                          + term4ik*rr7ik + term5ik*rr9ik;
                    }

                    // find standard multipole intermediates and energy value
                    else {
                        term1 = ci*ck;
                        term2 = ck*dir - ci*dkr + dik;
                        term3 = ci*qkr + ck*qir - dir*dkr + 2.*(dkqi-diqk+qiqk);
                        term4 = dir*qkr - dkr*qir - 4.*qik;
                        term5 = qir*qkr;
                        e = term1*rr1 + term2*rr3 + term3*rr5 + term4*rr7 + term5*rr9;
                    }

                    // increment the overall multipole energy components
                    if (use_group) e = e * fgrp;
                    if (e != 0.) {
                        nem++;
                        em += e;
                        aem[i] = aem[i] + 0.5*e;
                        aem[k] = aem[k] + 0.5*e;
                        if (molcule[i] != molcule[k]) {
                            einter += e;
                        }
                    }

                    // print message if the energy of this interaction is large
                    huge = (std::abs(e) > 100.);
                    if ((debug and e!=0.) or (verbose and huge)) {
                        if (header) {
                            header = false;
                            printf("\n Individual Atomic Multipole Interactions :\n\n");
                            printf(" Type              Atom Names               Distance        Energy\n\n");
                        }
                        printf(" Mpole     %7d-%-3s%7d-%-3s         %10.4f  %12.4f\n", i+1,name[i].c_str(),k+1,name[k].c_str(),r,e);
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

void empole_b(calcMode calculationMode){}
void empole_c(calcMode calculationMode){}
void empole_d(calcMode calculationMode){}
}
