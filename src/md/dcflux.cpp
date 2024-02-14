// Author: Moses KJ Chung
// Year:   2024

#include "angbnd.h"
#include "atoms.h"
#include "bndstr.h"
#include "bound.h"
#include "cflux.h"
#include "dcflux.h"
#include "libfunc.h"
#include "sizes.h"

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  dcflux  --  charge flux gradient chain rule  //
//                                               //
///////////////////////////////////////////////////

// "dcflux" takes as input the electrostatic potential at each
// atomic site and calculates gradient chain rule terms due to
// charge flux coupling with bond stretching and angle bending
//
// literature reference:
//
// C. Liu, J.-P. Piquemal and P. Ren, "Implementation of Geometry-
// Dependent Charge Flux into the Polarizable AMOEBA+ Potential",
// Journal of Physical Chemistry Letters, 11, 419-426 (2020)

template <CalcMode CalculationMode>
void dcflux(const real* cfpot, real* de, real (&virial)[3][3])
{
    int ia,ib,ic;
    real xa,ya,za;
    real xb,yb,zb;
    real xc,yc,zc;
    real xi,yi,zi;
    real xab,yab,zab;
    real xcb,ycb,zcb;
    real rab,rab2,rab3;
    real rcb,rcb2,rcb3;
    real dpot,dpota,dpotc;
    real fx,fy,fz;
    real fxa1,fya1,fza1;
    real fxb1,fyb1,fzb1;
    real fxc1,fyc1,fzc1;
    real fxa2,fya2,fza2;
    real fxb2,fyb2,fzb2;
    real fxc2,fyc2,fzc2;
    real pb,pb1,pb2;
    real pa1,pa2;
    real eps,dot;
    real rabc,dra3,drc3;
    real term,fterm;
    real termxa,termxc;
    real termya,termyc;
    real termza,termzc;

    // choose calculation mode
    constexpr CalcFlag flags = getCalculationFlags<CalculationMode>();
    constexpr bool do_v = flags.do_virial;

    // set pointers for OpenMP
    real* decfxP = decfx.ptr();
    real* decfyP = decfy.ptr();
    real* decfzP = decfz.ptr();

    // zero out the charge flux correction forces
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        decfxP[i] = 0.;
        decfyP[i] = 0.;
        decfzP[i] = 0.;
    }

    // calculate the charge flux forces due to bond stretches
    #pragma omp parallel for default(private)   \
    shared(nbond,x,y,z,ibnd,bflx,cfpot)         \
    reduction(+:decfxP[:n],decfyP[:n],decfzP[:n])
    for (int i = 0; i < nbond; i++) {
        ia = ibnd[i][0];
        ib = ibnd[i][1];
        pb = bflx[i];
        xa = x[ia];
        ya = y[ia];
        za = z[ia];
        xb = x[ib];
        yb = y[ib];
        zb = z[ib];
        xab = xa - xb;
        yab = ya - yb;
        zab = za - zb;
        rab2 = xab*xab + yab*yab + zab*zab;
        dpot = cfpot[ib] - cfpot[ia];
        pb = pb * dpot / REAL_SQRT(rab2);
        fx = pb * xab;
        fy = pb * yab;
        fz = pb * zab;
        decfxP[ia] += fx;
        decfyP[ia] += fy;
        decfzP[ia] += fz;
        decfxP[ib] -= fx;
        decfyP[ib] -= fy;
        decfzP[ib] -= fz;
    }

    // calculate the charge flux forces due to angle bends
    #pragma omp parallel for default(private)    \
    shared(nangle,x,y,z,iang,aflx,abflx,cfpot)   \
    reduction(+:decfxP[:n],decfyP[:n],decfzP[:n])
    for (int i = 0; i < nangle; i++) {
        ia = iang[i][0];
        ib = iang[i][1];
        ic = iang[i][2];
        pa1 = aflx[i][0];
        pa2 = aflx[i][1];
        pb1 = abflx[i][0];
        pb2 = abflx[i][1];
        xa = x[ia];
        ya = y[ia];
        za = z[ia];
        xb = x[ib];
        yb = y[ib];
        zb = z[ib];
        xc = x[ic];
        yc = y[ic];
        zc = z[ic];
        xab = xa - xb;
        yab = ya - yb;
        zab = za - zb;
        xcb = xc - xb;
        ycb = yc - yb;
        zcb = zc - zb;
        rab2 = xab*xab + yab*yab + zab*zab;
        rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
        rab  = REAL_SQRT(rab2);
        rcb  = REAL_SQRT(rcb2);

        if (rab>0 and rcb>0) {
            // get terms corresponding to asymmetric bond stretches
            dpota = cfpot[ia] - cfpot[ib];
            dpotc = cfpot[ic] - cfpot[ib];
            pb1 = dpota * pb1;
            pb2 = dpotc * pb2;
            fxa1 = pb2 * xab/rab;
            fya1 = pb2 * yab/rab;
            fza1 = pb2 * zab/rab;
            fxc1 = pb1 * xcb/rcb;
            fyc1 = pb1 * ycb/rcb;
            fzc1 = pb1 * zcb/rcb;
            fxb1 = -fxa1 - fxc1;
            fyb1 = -fya1 - fyc1;
            fzb1 = -fza1 - fzc1;

            // get terms corresponding to bond angle bending
            rabc = rab * rcb;
            rab3 = rab2 * rab;
            rcb3 = rcb2 * rcb;
            dot = xab*xcb + yab*ycb + zab*zcb;
            dra3 = dot / (rab3*rcb);
            drc3 = dot / (rab*rcb3);
            term = -rabc / REAL_SQRT(rab2*rcb2-dot*dot);
            fterm = term * (dpota*pa1+dpotc*pa2);
            termxa = xcb/rabc - xab*dra3;
            termya = ycb/rabc - yab*dra3;
            termza = zcb/rabc - zab*dra3;
            termxc = xab/rabc - xcb*drc3;
            termyc = yab/rabc - ycb*drc3;
            termzc = zab/rabc - zcb*drc3;
            fxa2 = fterm * termxa;
            fya2 = fterm * termya;
            fza2 = fterm * termza;
            fxc2 = fterm * termxc;
            fyc2 = fterm * termyc;
            fzc2 = fterm * termzc;
            fxb2 = -fxa2 - fxc2;
            fyb2 = -fya2 - fyc2;
            fzb2 = -fza2 - fzc2;
            decfxP[ia] += fxa1 + fxa2;
            decfyP[ia] += fya1 + fya2;
            decfzP[ia] += fza1 + fza2;
            decfxP[ib] += fxb1 + fxb2;
            decfyP[ib] += fyb1 + fyb2;
            decfzP[ib] += fzb1 + fzb2;
            decfxP[ic] += fxc1 + fxc2;
            decfyP[ic] += fyc1 + fyc2;
            decfzP[ic] += fzc1 + fzc2;
        }
    }

    // modify the gradient and virial for charge flux
    #pragma omp parallel for default(private)   \
    shared(n,x,y,z,decfx,decfy,decfz)           \
    reduction(+:virial,de[:3*n])
    for (int i = 0; i < n; i++) {
        fx = decfx[i];
        fy = decfy[i];
        fz = decfz[i];
        de[3*i]   += fx;
        de[3*i+1] += fy;
        de[3*i+2] += fz;

        if constexpr(do_v) {
            xi = x[i];
            yi = y[i];
            zi = z[i];
            real vxx = xi * fx;
            real vxy = yi * fx;
            real vxz = zi * fx;
            real vyy = yi * fy;
            real vyz = zi * fy;
            real vzz = zi * fz;
            virial[0][0] += vxx;
            virial[0][1] += vxy;
            virial[0][2] += vxz;
            virial[1][0] += vxy;
            virial[1][1] += vyy;
            virial[1][2] += vyz;
            virial[2][0] += vxz;
            virial[2][1] += vyz;
            virial[2][2] += vzz;
        }
    }
}

// explicit instatiation
template void dcflux<CalcMode::Gradient>(const real* cfpot, real* de, real (&virial)[3][3]);
template void dcflux<CalcMode::Virial>(const real* cfpot, real* de, real (&virial)[3][3]);
}
