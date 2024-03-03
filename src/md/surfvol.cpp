// Author: Moses KJ Chung
// Year:   2024

#include "alphmol.h"
#include "atoms.h"
#include "fatal.h"
#include "libfunc.h"
#include "mathConst.h"
#include "neigh.h"
#include "surfvol.h"

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  surfvol  --  compute surface area and volume  //
//                                                //
////////////////////////////////////////////////////

// "surfvol" uses Richmond's method to compute the
// surface area and volume
//
// literature references:
//
// T. J. Richmond, "Solvent Accessible Surface Area and
// Excluded Volume in Proteins", Journal of Molecular Biology,
// 178, 63-89 (1984)

inline real gaussbonnets(int jb, real arclen, real exang, real ri2);

void surfvol(bool deriv)
{
    int ncount;
    int jb,kb;
    real xi,yi,zi,ri,ri2;
    real xr,yr,zr,rk,rk2;
    real r,r2,rr,rr2;
    real s,v;
    real cos_th;
    real arclen,arcsum,exang;
    
    maxslst = 500;

    constexpr real pix2 = 2 * pi;
    constexpr real pix4 = 4 * pi;
    constexpr real pix4d3 = 4 * pi / 3;
    constexpr real pid2 = 0.5 * pi;
    constexpr real delta = 1e-8;
    constexpr real delta2 = delta * delta;
    constexpr real eps = 1e-8;

    // initialize surface area and volume
    tsurf = 0;
    tvol = 0;
    for (int i = 0; i < n; i++) {
        surf[i] = 0;
        vol[i] = 0;
    }
    if (deriv) {
        for (int i = 0; i < 3*n; i++) {
            dsurf[i] = 0;
            dvol[i] = 0;
        }
    }

    // allocate neighbor list
    nslst.allocate(n);
    slst.allocate(n*maxslst);

    // initialize neighbor list
    for (int i = 0; i < n; i++) {
        nslst[i] = 0;
    }

    // construct neighbor list
    for (int i = 0; i < n; i++) {
        xi = x[i];
        yi = y[i];
        zi = z[i];
        ri = radii[i];
        if (ri == 0.) continue;
        ncount = 0;
        bool buried = false;

        // evaluate all other spheres
        for (int k = 0; k < n; k++) {
            if (i == k) continue;
            xr = x[k] - xi;
            yr = y[k] - yi;
            zr = z[k] - zi;
            rk = radii[k];
            if (rk == 0.) continue;
            r2 = xr*xr + yr*yr + zr*zr;
            r = REAL_SQRT(r2);
            rr = ri + rk;
            rr2 = rr * rr;

            // check if i sphere is buried in k sphere
            if (rk > (r+ri)) {
                buried = true;
                break;
            }

            // include in neighborlist if (ri + rk)**2 is
            // greater than x**2 + y**2 + z**2
            if (rr2 > r2) {
                slst[i*maxslst+ncount] = k;
                ncount++;
            }
        }

        // check if there are too many overlapping spheres
        if (ncount > maxslst) {
            printf("\n SURFVOL  --  Too many intersections; Increase MAXSLST\n");
            fatal();
        }

        // if buried, set ncount = -1
        if (buried) ncount = -1;

        // store the number of neighbors
        nslst[i] = ncount;
    }

    // compute surface area and volume
    for (int i = 0; i < n; i++) {

        // buried atom case
        ncount = nslst[i];
        if (ncount == -1) continue;

        xi = x[i];
        yi = y[i];
        zi = z[i];
        ri = radii[i];
        ri2 = ri * ri;
        ncount = nslst[i];
        
        // surface area and volume of sphere
        s = pix4 * ri2;
        v = pix4d3 * ri2*ri;

        // if no neighbors
        if (ncount == 0) {
            surf[i] = s;
            vol[i] = v;
            tsurf += s;
            tvol += v;
        }

        // initialize arc length and external angle
        arclen = 0;
        exang = 0;

        // 
        jb = 0;
        kb = 0;

        // if one neighbor
        if (ncount == 1) {
            int k = slst[i*maxslst];
            xr = x[k] - xi;
            yr = y[k] - yi;
            zr = z[k] - zi;
            rk = radii[k];
            rk2 = rk * rk;
            r2 = xr*xr + yr*yr + zr*zr;
            r = REAL_SQRT(r2);

            // if i sphere buries k sphere
            if (ri > (r+rk)) {
                surf[i] = s;
                vol[i] = v;
                tsurf += s;
                tvol += v;
                continue;
            }

            // compute intersection angle theta
            cos_th = (r2 + ri2 - rk2) / (2 * r * ri);
            cos_th = REAL_MIN(1, REAL_MAX(-1, cos_th));
            
            arcsum = pix2;
            arclen += cos_th*arcsum;
            jb++;
            
            // compute surface area and volume using Gauss-Bonnet Theorem
            s = gaussbonnets(jb, arclen, exang, ri2);

            // increment surface area and volume
            surf[i] = s;
            tsurf += s;
        }
    }
}


inline real gaussbonnets(int jb, real arclen, real exang, real ri2)
{
    constexpr real pix2 = 2 * pi;
    constexpr real pix4 = 4 * pi;

    real s = jb*pix2 + exang + arclen;
    s = fmod(s, pix4);
    s *= ri2;
    return s;
}
}
