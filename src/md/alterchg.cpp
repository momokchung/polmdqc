// Author: Moses KJ Chung
// Year:   2024

#include "alterchg.h"
#include "angbnd.h"
#include "atmlst.h"
#include "atoms.h"
#include "bndstr.h"
#include "bound.h"
#include "cflux.h"
#include "chgpen.h"
#include "inform.h"
#include "libfunc.h"
#include "mathConst.h"
#include "mplpot.h"
#include "mpole.h"
#include "sizes.h"

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  alterchg  --  modification of partial charges  //
//                                                 //
/////////////////////////////////////////////////////

// "alterchg" calculates the change in atomic partial charge or
// monopole values due to bond and angle charge flux coupling

// literature reference:

// C. Liu, J.-P. Piquemal and P. Ren, "Implementation of Geometry-
// Dependent Charge Flux into the Polarizable AMOEBA+ Potential",
// Journal of Physical Chemistry Letters, 11, 419-426 (2020)

void alterchg()
{
    // zero out the change in charge value at each site
    for (int i = 0; i < n; i++) {
        pdelta[i] = 0.;
    }

    // find charge modifications due to charge flux
    bndchg(pdelta.ptr());
    angchg(pdelta.ptr());

    // alter monopoles and charge penetration for charge flux
    #pragma omp parallel for default(private)   \
    shared(n,pdelta,pole,mono0,pval,pval0)
    for (int i = 0; i < n; i++) {
        pole[i][0] = mono0[i] + pdelta[i];
        pval[i] = pval0[i] + pdelta[i];
    }
}

/////////////////////////////////////////////////////
//                                                 //
//  bndchg  --  charge flux bond stretch coupling  //
//                                                 //
/////////////////////////////////////////////////////

// "bndchg" computes modifications to atomic partial charges or
// monopoles due to bond stretch using a charge flux formulation

void bndchg(real* pdelta)
{
    int ia,ib;
    real xab,yab,zab;
    real rab,rab0;
    real pb,dq;

    // loop over all the bond distances in the system
    #pragma omp parallel for default(private)   \
    shared(nbond,x,y,z,ibnd,bl,bflx)            \
    reduction(+:pdelta[:n])
    for (int i = 0; i < nbond; i++) {
        ia = ibnd[i][0];
        ib = ibnd[i][1];
        pb = bflx[i];

        // compute the bond length value for the current bond
        xab = x[ia] - x[ib];
        yab = y[ia] - y[ib];
        zab = z[ia] - z[ib];
        rab = REAL_SQRT(xab*xab + yab*yab + zab*zab);

        // find the charge flux increment for the current bond
        rab0 = bl[i];
        dq = pb * (rab-rab0);
        pdelta[ia] -= dq;
        pdelta[ib] += dq;
    }
}

///////////////////////////////////////////////////
//                                               //
//  angchg  --  charge flux angle bend coupling  //
//                                               //
///////////////////////////////////////////////////

// "angchg" computes modifications to atomic partial charges or
// monopoles due to angle bending using a charge flux formulation

void angchg(real* pdelta)
{
    int ia,ib,ic;
    real angle;
    real rab,rcb;
    real xia,yia,zia;
    real xib,yib,zib;
    real xic,yic,zic;
    real xab,yab,zab;
    real xcb,ycb,zcb;
    real dot,cosine;
    real pa1,pa2;
    real pb1,pb2;
    real theta0;
    real rab0,rcb0;
    real dq1,dq2;

    // loop over all the bond angles in the system
    #pragma omp parallel for default(private)             \
    shared(nangle,x,y,z,iang,aflx,abflx,anat,bl,balist)   \
    reduction(+:pdelta[:n])
    for (int i = 0; i < nangle; i++) {
        ia = iang[i][0];
        ib = iang[i][1];
        ic = iang[i][2];
        pa1 = aflx[i][0];
        pa2 = aflx[i][1];
        pb1 = abflx[i][0];
        pb2 = abflx[i][1];

        // calculate the angle values and included bond lengths
        xia = x[ia];
        yia = y[ia];
        zia = z[ia];
        xib = x[ib];
        yib = y[ib];
        zib = z[ib];
        xic = x[ic];
        yic = y[ic];
        zic = z[ic];
        xab = xia - xib;
        yab = yia - yib;
        zab = zia - zib;
        xcb = xic - xib;
        ycb = yic - yib;
        zcb = zic - zib;
        rab = REAL_SQRT(xab*xab+yab*yab+zab*zab);
        rcb = REAL_SQRT(xcb*xcb+ycb*ycb+zcb*zcb);
        if (rab>0 and rcb>0) {
            dot = xab*xcb + yab*ycb + zab*zcb;
            cosine = dot / (rab*rcb);
            cosine = REAL_MIN((real)1,REAL_MAX((real)-1,cosine));
            angle = radian * REAL_ACOS(cosine);

            // find the charge flux increment for the current angle
            theta0 = anat[i];
            rab0 = bl[balist[i][0]];
            rcb0 = bl[balist[i][1]];
            dq1 = pb1*(rcb-rcb0) + pa1*(angle-theta0)/radian;
            dq2 = pb2*(rab-rab0) + pa2*(angle-theta0)/radian;
            pdelta[ia] += dq1;
            pdelta[ib] -= dq1 + dq2;
            pdelta[ic] += dq2;
        }
    }
}
}
